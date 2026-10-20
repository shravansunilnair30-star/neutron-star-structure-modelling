"""Convergence study for the CompOSE tabulated-EOS TOV solver.
Requires data/eos.thermo and data/eos.nb downloaded from CompOSE RG(SLY4).
"""
from pathlib import Path
import csv
import numpy as np
from neutron_star_compose import read_compose_table, build_eos, integrate_star, interp_radius_at_mass

ROOT = Path(__file__).resolve().parent
DATA = ROOT/'data'
OUT = ROOT/'outputs'

def run_case(pcs, p_all, interp, rtol, max_step):
    rows = [r for pc in pcs if (r := integrate_star(pc, p_all, interp, rtol, max_step)) is not None]
    if len(rows) < 20:
        raise RuntimeError(f'Only {len(rows)} integrations succeeded for rtol={rtol}, max_step={max_step}.')
    masses = np.array([r['mass_solar'] for r in rows]); radii = np.array([r['radius_km'] for r in rows])
    imax = int(np.argmax(masses))
    return {'rtol':rtol,'max_step_cm':max_step,'successful_models':len(rows),
            'max_sampled_mass_solar':float(masses[imax]),'radius_at_max_mass_km':float(radii[imax]),
            'radius_at_1p4_solar_mass_km':float(interp_radius_at_mass(rows))}

def main():
    OUT.mkdir(exist_ok=True)
    ptab, etab = read_compose_table(DATA/'eos.thermo', DATA/'eos.nb')
    p_all, eps_all, interp, p_causal_max, cs2, causal = build_eos(ptab, etab)
    pcs = np.geomspace(max(ptab[0]*1.01,1e28), min(p_causal_max*.999,ptab[-1]*.999), 140)
    configs = [(2e-7,2e4),(2e-8,1e4),(2e-9,5e3),(2e-10,2.5e3)]
    results = [run_case(pcs,p_all,interp,rtol,step) for rtol,step in configs]
    base = results[1]
    for row in results:
        for key in ('max_sampled_mass_solar','radius_at_max_mass_km','radius_at_1p4_solar_mass_km'):
            row[key+'_relative_change_vs_reference_pct'] = (row[key]-base[key])/base[key]*100
    with (OUT/'compose_sly4_convergence.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=results[0].keys());w.writeheader();w.writerows(results)
    summary=['CompOSE RG(SLY4) tabulated-EOS convergence study',
             'Reference configuration: rtol=2e-8, max_step=1e4 cm; same central-pressure grid in every run.']
    for r in results:
        summary.append(f"rtol={r['rtol']:.1e}, max_step={r['max_step_cm']:.0f} cm: Mmax={r['max_sampled_mass_solar']:.7f} Msun, Rmax={r['radius_at_max_mass_km']:.6f} km, R1.4={r['radius_at_1p4_solar_mass_km']:.6f} km; ΔM={r['max_sampled_mass_solar_relative_change_vs_reference_pct']:+.5f}%, ΔRmax={r['radius_at_max_mass_km_relative_change_vs_reference_pct']:+.5f}%, ΔR1.4={r['radius_at_1p4_solar_mass_km_relative_change_vs_reference_pct']:+.5f}%")
    (OUT/'compose_sly4_convergence_summary.txt').write_text('\n'.join(summary)+'\n')
    print('\n'.join(summary))
if __name__=='__main__': main()
