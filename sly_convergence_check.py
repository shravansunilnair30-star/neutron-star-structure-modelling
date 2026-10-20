"""Numerical sensitivity check for the SLy piecewise-polytropic TOV solver."""
import csv
from pathlib import Path
import numpy as np
from neutron_star_sly import integrate_star

OUT=Path(__file__).resolve().parent/'outputs'; OUT.mkdir(exist_ok=True)
pcs=np.logspace(34.0,36.5,70)
configs=[('baseline',2e-8,2e4),('tight_tolerance',2e-9,2e4),('half_max_step',2e-8,1e4),('tight_and_half_step',2e-9,1e4)]
rows=[]
for name,rtol,max_step in configs:
    vals=[r for pc in pcs if (r:=integrate_star(pc,rtol=rtol,max_step=max_step)) is not None]
    i=int(np.argmax([r['mass_solar'] for r in vals])); peak=vals[i]
    rows.append({'configuration':name,'rtol':rtol,'max_step_cm':max_step,'successful_models':len(vals),'max_sampled_mass_solar':peak['mass_solar'],'radius_at_max_km':peak['radius_km']})
with (OUT/'sly_convergence_test.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
base=rows[0]
lines=['SLy piecewise-polytrope numerical sensitivity check','Comparisons are relative to baseline; this tests solver stability, not physical EOS accuracy.','']
for r in rows:
    dm=100*(r['max_sampled_mass_solar']/base['max_sampled_mass_solar']-1)
    dr=100*(r['radius_at_max_km']/base['radius_at_max_km']-1)
    lines.append(f"{r['configuration']}: Mmax(sampled)={r['max_sampled_mass_solar']:.7f} Msun, R={r['radius_at_max_km']:.5f} km, ΔM={dm:+.6g}%, ΔR={dr:+.6g}%")
(OUT/'sly_convergence_test_summary.txt').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
