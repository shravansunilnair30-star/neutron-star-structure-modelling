"""Physical checks for the Read et al. SLy piecewise-polytropic TOV approximation.

This validates the *approximation* against published CompOSE RG(SLY4) summary
properties. It does NOT claim that the CompOSE tabulated EOS was loaded by this
script. All outputs are generated locally from neutron_star_sly.py.
"""
from __future__ import annotations
import csv
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import neutron_star_sly as eos

OUT = Path(__file__).resolve().parent / 'outputs'
OUT.mkdir(exist_ok=True)
# Rounded reference properties published on the official CompOSE RG(SLY4) page.
REF = {'Mmax_Msun': 2.06, 'R_at_Mmax_km': 10.02, 'R14_km': 11.70}


def cs2_over_c2(rho: float, segment: int) -> float:
    """dP/dε for ε=(1+A)rho*c^2 + P/(Gamma-1), in one polytropic segment."""
    g = eos.GAMMAS[segment]
    k = eos.K[segment]
    dp_drho = g * k * rho ** (g - 1)
    de_drho = (1 + eos.A[segment]) * eos.C**2 + g/(g-1) * k * rho**(g-1)
    return float(dp_drho / de_drho)


def segment_for_density(rho: float) -> int:
    if rho < eos.RHO0: return 0
    if rho < eos.RHO1: return 1
    if rho < eos.RHO2: return 2
    return 3


def max_interior_cs2(central_rho: float) -> tuple[float, float, int]:
    """Maximum c_s^2/c^2 over every EOS segment sampled by a star interior.

    The analytic c_s^2 expression is monotone increasing within each segment;
    evaluate both ends of every segment intersecting [surface, central density].
    Returns maximum, density at maximum, and segment index.
    """
    bounds = [1e10, eos.RHO0, eos.RHO1, eos.RHO2, central_rho]
    best = (-np.inf, np.nan, -1)
    for i in range(4):
        lo = bounds[i]
        hi = min(bounds[i+1], central_rho)
        if hi <= lo:
            continue
        for rho in (lo, hi):
            value = cs2_over_c2(rho, i)
            if value > best[0]: best = (value, rho, i)
    # Include all boundaries from the lower-density side too, to capture derivative jumps.
    for i, rho in enumerate([eos.RHO0, eos.RHO1, eos.RHO2]):
        if rho <= central_rho:
            for seg in (i, i+1):
                value = cs2_over_c2(rho, seg)
                if value > best[0]: best = (value, rho, seg)
    return float(best[0]), float(best[1]), int(best[2])


def interpolate_r14(rows):
    rows = sorted(rows, key=lambda r: r['central_pressure_dyn_cm2'])
    # First crossing on the low-central-pressure branch.
    for a, b in zip(rows[:-1], rows[1:]):
        if a['mass_solar'] <= 1.4 <= b['mass_solar'] and b['mass_solar'] > a['mass_solar']:
            f = (1.4-a['mass_solar'])/(b['mass_solar']-a['mass_solar'])
            return a['radius_km'] + f*(b['radius_km']-a['radius_km'])
    return float('nan')


def run_configuration(rtol, max_step):
    # Find where the analytic EOS first reaches c_s = c in its high-density segment.
    def causal_eq(rho): return cs2_over_c2(rho, 3) - 1.0
    rho_causal = brentq(causal_eq, eos.RHO2, 4e15, xtol=1e-3)
    p_causal = eos.K[3] * rho_causal**eos.GAMMAS[3]
    # Keep a tiny safety margin below the causal boundary to avoid rounding over it.
    pc = p_causal * (1 - 1e-7)
    star = eos.integrate_star(pc, rtol=rtol, max_step=max_step)
    if star is None: raise RuntimeError('TOV integration failed at causal boundary.')
    max_cs2, rho_at_max, seg = max_interior_cs2(star['central_density_g_cm3'])
    # Dense low-branch sequence for R(1.4) and sampled causal maximum.
    pcs = np.geomspace(1e34, pc, 240)
    rows = [s for p in pcs if (s := eos.integrate_star(p, rtol=rtol, max_step=max_step)) is not None]
    imax = int(np.argmax([r['mass_solar'] for r in rows]))
    peak = rows[imax]
    return {
        'rtol': rtol, 'max_step_cm': max_step, 'rho_causal_g_cm3': rho_causal,
        'p_causal_dyn_cm2': p_causal, 'causal_boundary_mass_solar': star['mass_solar'],
        'causal_boundary_radius_km': star['radius_km'], 'max_interior_cs2_over_c2': max_cs2,
        'density_at_max_cs2_g_cm3': rho_at_max, 'segment_at_max_cs2': seg,
        'sampled_causal_max_mass_solar': peak['mass_solar'], 'sampled_causal_max_radius_km': peak['radius_km'],
        'r14_km': interpolate_r14(rows), 'successful_models': len(rows), 'rows': rows,
    }


def main():
    configs = [('baseline', 1e-8, 2e4), ('tight_tol', 1e-9, 2e4),
               ('half_step', 1e-8, 1e4), ('tight_and_half_step', 1e-9, 1e4)]
    results = []
    for name, rtol, step in configs:
        r = run_configuration(rtol, step); r['configuration'] = name; results.append(r)
    base = results[0]
    # Save all causal-branch sequence results for the baseline configuration.
    with (OUT/'physical_validation_causal_mass_radius.csv').open('w', newline='') as f:
        rows = base.pop('rows')
        w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    # Restore other rows key removal irrelevant; baseline has its values already saved.
    fields = [k for k in results[0].keys() if k != 'rows']
    with (OUT/'physical_validation_convergence.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for r in results: w.writerow({k:r.get(k, '') for k in fields})
    # Main values based on baseline and official rounded benchmarks.
    m, rmax, r14 = base['causal_boundary_mass_solar'], base['causal_boundary_radius_km'], base['r14_km']
    errors = {'Mmax': 100*(m/REF['Mmax_Msun']-1),
              'R_at_Mmax': 100*(rmax/REF['R_at_Mmax_km']-1),
              'R14': 100*(r14/REF['R14_km']-1)}
    max_conv_m = max(abs(100*(r['causal_boundary_mass_solar']/m-1)) for r in results)
    max_conv_r = max(abs(100*(r['causal_boundary_radius_km']/rmax-1)) for r in results)
    causal_pass = all(r['max_interior_cs2_over_c2'] <= 1 + 1e-8 for r in results)
    convergence_pass = max(max_conv_m, max_conv_r) < 0.1
    reference_pass = max(abs(v) for v in errors.values()) < 2.0
    lines = [
        'PHYSICAL VALIDATION RESULTS — SLy PIECEWISE-POLYTROPE APPROXIMATION',
        '',
        'Model identity: Read et al. (2009) SLy piecewise-polytropic approximation with matched low-density polytropic crust.',
        'Important: this run does NOT load the CompOSE eos.thermo table. It compares this analytic approximation with the rounded reference properties on the official CompOSE RG(SLY4) page.',
        '',
        'EOS causality test (full interior):',
        f'  Density where analytic EOS reaches c_s=c: {base["rho_causal_g_cm3"]:.6e} g/cm^3',
        f'  Corresponding central pressure: {base["p_causal_dyn_cm2"]:.6e} dyn/cm^2',
        f'  Maximum c_s^2/c^2 over causal-boundary stellar interior: {base["max_interior_cs2_over_c2"]:.10f}',
        f'  Density of maximum c_s^2/c^2: {base["density_at_max_cs2_g_cm3"]:.6e} g/cm^3; EOS segment index: {base["segment_at_max_cs2"]}',
        f'  Causality check: {"PASS" if causal_pass else "FAIL"}',
        '',
        'Results at causal boundary (tiny 1e-7 fractional safety margin):',
        f'  Mass: {m:.7f} M_sun',
        f'  Radius: {rmax:.6f} km',
        f'  R(1.4 M_sun), interpolated from causal sequence: {r14:.6f} km',
        '',
        'Official CompOSE RG(SLY4) page reference (rounded):',
        f'  Mmax = {REF["Mmax_Msun"]:.3f} M_sun; R at Mmax = {REF["R_at_Mmax_km"]:.3f} km; R(1.4) = {REF["R14_km"]:.3f} km',
        f'  Relative differences: Mmax {errors["Mmax"]:+.3f}%; R at Mmax {errors["R_at_Mmax"]:+.3f}%; R(1.4) {errors["R14"]:+.3f}%',
        f'  Scalar-reference comparison (<2% all metrics): {"PASS" if reference_pass else "FAIL"}',
        '',
        'Numerical convergence:',
        f'  Largest absolute change across tested tolerances/step sizes: mass {max_conv_m:.6g}%; radius {max_conv_r:.6g}%',
        f'  Convergence check (<0.1%): {"PASS" if convergence_pass else "FAIL"}',
        '',
        'OVERALL INTERPRETATION:',
        f'  Approximation-level checks: {"PASS" if (causal_pass and convergence_pass and reference_pass) else "PARTIAL / FAIL"}.',
        '  This is encouraging evidence that the causal branch of the analytic SLy approximation is numerically stable and close to the published RG(SLY4) global properties.',
        '  It is NOT a full validation of a tabulated SLy4 EOS implementation: the original CompOSE table was not loaded, and the compared model uses an analytic piecewise-polytropic fit with an approximate crust. Do not describe it as an exact reproduction of tabulated SLy4.',
        '  A full tabulated-EOS validation still requires running neutron_star_compose.py with the official cold CompOSE files (eos.thermo and eos.nb), checking its full EOS causality curve and convergence, and comparing its computed sequence to the official eos.mr curve.',
    ]
    text='\n'.join(lines)+'\n'
    (OUT/'physical_validation_summary.txt').write_text(text)
    print(text)

if __name__ == '__main__': main()
