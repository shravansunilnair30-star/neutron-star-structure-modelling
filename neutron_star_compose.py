"""TOV solver driven by an official CompOSE cold beta-equilibrium EOS table.

Expected inputs (download RG(SLY4) from https://compose.obspm.fr/eos/134):
  data/eos.thermo and data/eos.nb

The script reads the CompOSE standard eos.thermo columns, constructs P(epsilon),
checks thermodynamic monotonicity and causality, integrates TOV stars only up to
the causal pressure limit, and compares global properties with the reference
values listed on the CompOSE RG(SLY4) page. It does not silently treat a failed
check as a pass.

CompOSE eos.thermo conventions are documented in the CompOSE reference manual.
"""
from __future__ import annotations
import argparse, csv, math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp

G = 6.67430e-8
C = 2.99792458e10
MSUN = 1.98847e33
MEV_FM3_TO_CGS = 1.602176634e33
KM = 1.0e5
# Approximate reference values listed on CompOSE RG(SLY4) model page.
REF_MMAX = 2.06
REF_RMAX_KM = 10.02
REF_R14_KM = 11.70


def read_compose_table(thermo_path: Path, nb_path: Path):
    """Return sorted pressure and energy density arrays in cgs.

    CompOSE eos.thermo rows after the first header row are:
    iT, inb, iYq, P/nb, s/nb, mub/mn-1, muq/mn, mul/mn,
    f/(nb mn)-1, e/(nb mn)-1, Nadd, [additional quantities].
    The first three entries in eos.thermo's first line are mn, mp (MeV),
    and a leptons flag. eos.nb gives the baryon-density grid in fm^-3.
    """
    raw = np.loadtxt(thermo_path, skiprows=1)
    if raw.ndim == 1:
        raw = raw[None, :]
    header = np.loadtxt(thermo_path, max_rows=1)
    if len(header) < 2:
        raise ValueError('CompOSE eos.thermo first line must contain neutron/proton masses and lepton flag.')
    mn_mev = float(header[0])
    # eos.nb: first two lines are min/max indices, following lines are values.
    nb_lines = [x.strip() for x in nb_path.read_text().splitlines() if x.strip() and not x.lstrip().startswith('#')]
    if len(nb_lines) < 3:
        raise ValueError('eos.nb does not contain the expected index bounds and density grid.')
    nb_values = np.array([float(x.split()[0]) for x in nb_lines[2:]], dtype=float)
    nb_min_idx = int(float(nb_lines[0].split()[0]))
    # CompOSE index values are 1-based/explicit in eos.thermo. Map index to grid by offset.
    nb_idx = raw[:, 1].astype(int)
    offsets = nb_idx - nb_min_idx
    if offsets.min() < 0 or offsets.max() >= len(nb_values):
        # Some files use zero-based indices; try that convention only if consistent.
        offsets = nb_idx - 0
    if offsets.min() < 0 or offsets.max() >= len(nb_values):
        raise ValueError('Density indices in eos.thermo do not match eos.nb grid indices.')
    nb = nb_values[offsets]
    # Standard thermodynamic quantities are at columns 3..9; Q1=P/nb and Q7=e/(nb mn)-1.
    p_mev_fm3 = raw[:, 3] * nb
    eps_mev_fm3 = nb * mn_mev * (1.0 + raw[:, 9])
    # Cold beta-equilibrium tables should have one temperature and one charge-fraction index.
    keys = np.unique(raw[:, [0, 2]], axis=0)
    if len(keys) != 1:
        raise ValueError(f'Expected a 1D cold beta-equilibrium table, found {len(keys)} (T,Yq) index combinations.')
    order = np.argsort(p_mev_fm3)
    p = p_mev_fm3[order] * MEV_FM3_TO_CGS
    eps = eps_mev_fm3[order] * MEV_FM3_TO_CGS
    keep = np.isfinite(p) & np.isfinite(eps) & (p > 0) & (eps > 0)
    p, eps = p[keep], eps[keep]
    # Collapse any duplicate pressures to avoid invalid interpolation.
    p_unique, unique_idx = np.unique(p, return_index=True)
    eps = eps[unique_idx]
    if len(p_unique) < 30:
        raise ValueError(f'Only {len(p_unique)} unique valid EOS points were read; check table format.')
    if np.any(np.diff(p_unique) <= 0) or np.any(np.diff(eps) <= 0):
        raise ValueError('EOS pressure and energy density must be strictly increasing for this barotropic solver.')
    return p_unique, eps


def build_eos(p_table, eps_table):
    """Build log-log PCHIP epsilon(P), with a low-pressure power-law extension.

    The tabulated CompOSE EOS starts at finite density. The extension to P=1 dyn/cm²
    is fitted to the lowest 8 table points and is explicitly an extrapolation. Its
    effect should be checked by varying the number of points used in the fit.
    """
    nfit = min(8, len(p_table))
    slope, intercept = np.polyfit(np.log(p_table[:nfit]), np.log(eps_table[:nfit]), 1)
    # epsilon ~ P^slope; equivalent P ~ epsilon^gamma, gamma=1/slope.
    if not (0.0 < slope < 1.0):
        raise ValueError(f'Low-density extrapolation exponent {slope:.4g} is implausible; inspect EOS table.')
    p_floor = 1.0
    eps_floor = float(np.exp(intercept + slope*np.log(p_floor)))
    p_ext = np.geomspace(p_floor, p_table[0], 80, endpoint=False)
    eps_ext = np.exp(intercept + slope*np.log(p_ext))
    p_all = np.concatenate([p_ext, p_table])
    eps_all = np.concatenate([eps_ext, eps_table])
    # Deduplicate and interpolate log epsilon against log pressure.
    p_all, ix = np.unique(p_all, return_index=True)
    eps_all = eps_all[ix]
    interp = PchipInterpolator(np.log(p_all), np.log(eps_all), extrapolate=False)
    dlogeps_dlogp = interp.derivative()
    # c_s^2/c^2 = dP/depsilon = (P/epsilon)/(d ln epsilon/d ln P).
    cs2 = (p_all / eps_all) / dlogeps_dlogp(np.log(p_all))
    causal = np.isfinite(cs2) & (cs2 > 0) & (cs2 <= 1.0 + 1e-6)
    # Only consider the continuous causal interval from the low-pressure end.
    bad = np.flatnonzero(~causal)
    p_causal_max = float(p_all[bad[0]-1]) if len(bad) and bad[0] > 0 else (float(p_all[-1]) if not len(bad) else 0.0)
    return p_all, eps_all, interp, p_causal_max, cs2, causal


def epsilon_of_p(p, p_all, logeps_interp):
    p = float(np.clip(p, p_all[0], p_all[-1]))
    return float(np.exp(logeps_interp(np.log(p))))


def integrate_star(pc, p_all, logeps_interp, rtol=2e-8, max_step=1e4):
    p_surface = float(p_all[0])
    epsc = epsilon_of_p(pc, p_all, logeps_interp)
    r0 = 1.0
    m0 = 4*np.pi*r0**3*epsc/(3*C**2)
    # Central Taylor expansion for pressure at the starting radius.
    dpdr2 = -G*(epsc+pc)*(m0+4*np.pi*r0**3*pc/C**2)/(C**2*r0**2*max(1-2*G*m0/(r0*C**2), 1e-15))
    p0 = pc + 0.5*dpdr2*r0**2
    def rhs(r, y):
        m, p = y
        p = max(float(p), p_surface)
        eps = epsilon_of_p(p, p_all, logeps_interp)
        compact = max(1 - 2*G*m/(r*C**2), 1e-14)
        dm = 4*np.pi*r*r*eps/C**2
        dp = -G*(eps+p)*(m+4*np.pi*r**3*p/C**2)/(C**2*r*r*compact)
        return [dm, dp]
    def surface(r, y): return y[1]-p_surface
    surface.terminal = True
    surface.direction = -1
    sol = solve_ivp(rhs, (r0, 3e7), [m0, p0], method='DOP853', rtol=rtol,
                    atol=[1e15, max(p_surface*1e-5, 1e-4)], max_step=max_step, events=surface)
    if not sol.success or not len(sol.t_events[0]):
        return None
    R = float(sol.t_events[0][0]); M = float(sol.y_events[0][0][0])
    return {'central_pressure_dyn_cm2': float(pc), 'mass_solar': M/MSUN, 'radius_km': R/KM}


def interp_radius_at_mass(rows, target=1.4):
    # Find adjacent mass samples bracketing target on the low-central-pressure branch.
    rows = sorted(rows, key=lambda r: r['central_pressure_dyn_cm2'])
    for a, b in zip(rows[:-1], rows[1:]):
        if a['mass_solar'] <= target <= b['mass_solar'] and b['mass_solar'] > a['mass_solar']:
            f = (target-a['mass_solar'])/(b['mass_solar']-a['mass_solar'])
            return a['radius_km'] + f*(b['radius_km']-a['radius_km'])
    return float('nan')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-dir', type=Path, default=Path('data'), help='Directory containing CompOSE eos.thermo and eos.nb')
    ap.add_argument('--nmodels', type=int, default=180)
    ap.add_argument('--rtol', type=float, default=2e-8)
    ap.add_argument('--max-step-cm', type=float, default=1e4)
    args = ap.parse_args()
    thermo = args.data_dir/'eos.thermo'; nbfile = args.data_dir/'eos.nb'
    if not thermo.exists() or not nbfile.exists():
        raise SystemExit('Missing CompOSE input. Download the RG(SLY4) table from https://compose.obspm.fr/eos/134, extract eos.thermo and eos.nb into data/.')
    ptab, etab = read_compose_table(thermo, nbfile)
    p_all, eps_all, interp, p_causal_max, cs2, causal = build_eos(ptab, etab)
    out = Path(__file__).resolve().parent/'outputs'; out.mkdir(exist_ok=True)
    # Only integrate stars with central pressures below the first causality failure.
    if p_causal_max <= p_all[0]*10:
        raise RuntimeError('No meaningful causal EOS interval found. Check CompOSE columns and interpolation.')
    pcs = np.geomspace(max(ptab[0]*1.01, 1e28), min(p_causal_max*0.999, ptab[-1]*0.999), args.nmodels)
    rows = [x for pc in pcs if (x := integrate_star(pc, p_all, interp, args.rtol, args.max_step_cm)) is not None]
    if len(rows) < 20: raise RuntimeError(f'Only {len(rows)} successful stellar models.')
    with (out/'compose_sly4_mass_radius.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    masses = np.array([r['mass_solar'] for r in rows]); radii = np.array([r['radius_km'] for r in rows])
    imax = int(np.argmax(masses)); r14 = interp_radius_at_mass(rows)
    fig, ax = plt.subplots(figsize=(8, 5.5)); ax.plot(radii, masses, 'o-', ms=2.5, lw=1.3, label='CompOSE RG(SLY4) tabulated EOS')
    ax.scatter([radii[imax]], [masses[imax]], marker='*', s=130, label=f'Maximum causal sampled mass: {masses[imax]:.3f} M☉')
    if np.isfinite(r14): ax.axhline(1.4, ls='--', lw=1, label=f'R(1.4 M☉) = {r14:.3f} km')
    ax.set_xlabel('Radius (km)'); ax.set_ylabel('Gravitational mass (solar masses)'); ax.set_title('TOV mass–radius relation using tabulated CompOSE RG(SLY4) EOS'); ax.grid(True, alpha=.3); ax.legend(); fig.tight_layout()
    fig.savefig(out/'compose_sly4_mass_radius.png', dpi=220); plt.close(fig)
    # Table-wide EOS causality diagnostics; derivative computed from the PCHIP curve.
    d = interp.derivative()(np.log(p_all))
    cs2_all = (p_all/np.exp(interp(np.log(p_all)))) / d
    eos_rows = [{'pressure_dyn_cm2': float(p), 'energy_density_erg_cm3': float(e), 'cs2_over_c2': float(v), 'causal': bool(np.isfinite(v) and 0 < v <= 1+1e-6)} for p,e,v in zip(p_all, np.exp(interp(np.log(p_all))), cs2_all)]
    with (out/'compose_sly4_eos_causality.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=eos_rows[0].keys()); w.writeheader(); w.writerows(eos_rows)
    rel_m = (masses[imax]-REF_MMAX)/REF_MMAX*100
    rel_rmax = (radii[imax]-REF_RMAX_KM)/REF_RMAX_KM*100
    rel_r14 = (r14-REF_R14_KM)/REF_R14_KM*100 if np.isfinite(r14) else float('nan')
    verdict = 'PASS (within 2% of all three approximate CompOSE scalar benchmarks)' if max(abs(rel_m),abs(rel_rmax),abs(rel_r14)) <= 2 else 'PARTIAL / INVESTIGATE: at least one scalar benchmark differs by more than 2%'
    summary = f'''PHYSICAL VALIDATION: COMPOSE RG(SLY4) EOS\n\nInput files: {thermo.name}, {nbfile.name}\nEOS table points used: {len(ptab)}\nCausal pressure limit from interpolated EOS: {p_causal_max:.6e} dyn/cm^2\nEOS causality rows flagged: {int(np.sum(~np.array([x['causal'] for x in eos_rows])))} / {len(eos_rows)} (includes low-pressure extrapolation and high-density range)\nSolver: DOP853, rtol={args.rtol:g}, max_step={args.max_step_cm:g} cm\nSuccessful stellar models: {len(rows)} / {len(pcs)}\n\nRESULTS (central pressures restricted to causal EOS interval)\nMaximum sampled mass: {masses[imax]:.6f} M_sun at R={radii[imax]:.5f} km\nInterpolated R_1.4: {r14:.5f} km\n\nCOMPOSE PAGE REFERENCE VALUES (approximate, same named RG(SLY4) model)\nM_max = {REF_MMAX:.3f} M_sun; R_Mmax = {REF_RMAX_KM:.3f} km; R_1.4 = {REF_R14_KM:.3f} km\nRelative differences: M_max={rel_m:+.3f}%; R_Mmax={rel_rmax:+.3f}%; R_1.4={rel_r14:+.3f}%\n\nAUTOMATIC SCALAR-COMPARISON VERDICT: {verdict}\n\nLIMITATIONS / INTERPRETATION\n- The CompOSE page's published scalar values are rounded reference targets, not an independent reference mass-radius curve.\n- The low-pressure extension below the first tabulated EOS point is a fitted power law; repeat with different fit ranges and lower pressure floors to quantify surface-radius sensitivity.\n- Causality is checked from the interpolated EOS derivative, not by relying on central density alone. Inspect compose_sly4_eos_causality.csv.\n- This is a computational reproduction against a named EOS dataset, not proof that SLy4 is the true astrophysical EOS.\n- Do not claim full validation unless the causality file contains no unexamined violations over the density range used and the convergence study passes.\n'''
    (out/'compose_sly4_validation_summary.txt').write_text(summary)
    print(summary)

if __name__ == '__main__': main()
