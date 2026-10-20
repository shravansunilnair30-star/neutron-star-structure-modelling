"""TOV solver using the SLy piecewise-polytropic approximation of Read et al. (2009).
All quantities in cgs. This is a literature-calibrated analytic approximation, not the original tabulated SLy4 microphysical EOS.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

G=6.67430e-8; C=2.99792458e10; MSUN=1.98847e33; KM=1e5
# Read et al. (2009) SLy piecewise-polytrope parameters.
RHO0=10**14.3; RHO1=10**14.7; RHO2=1e15
P1=10**34.384
GAMMAS=[1.35692395,3.005,2.988,2.851]
# Segment 0 is an approximate crust matched to the first core segment at rho0.
K=[0.0]*4
K[1]=P1/RHO1**GAMMAS[1]
P0=K[1]*RHO0**GAMMAS[1]
K[0]=P0/RHO0**GAMMAS[0]
K[2]=K[1]*RHO1**(GAMMAS[1]-GAMMAS[2])
P2=K[2]*RHO2**GAMMAS[2]
K[3]=K[2]*RHO2**(GAMMAS[2]-GAMMAS[3])
# Additive constants in epsilon(rho) enforce continuity at each density boundary.
A=[0.0]*4
for i, rt in [(1,RHO0),(2,RHO1),(3,RHO2)]:
    prev=i-1
    p=K[i]*rt**GAMMAS[i]
    eps_prev=(1+A[prev])*rt*C**2+p/(GAMMAS[prev]-1)
    A[i]=eps_prev/(rt*C**2)-1-p/((GAMMAS[i]-1)*rt*C**2)

def rho_from_pressure(p):
    # Pressure thresholds at rho0, rho1 and rho2.
    if p < P0: i=0
    elif p < P1: i=1
    elif p < P2: i=2
    else: i=3
    return (p/K[i])**(1/GAMMAS[i]), i

def energy_density(p):
    rho,i=rho_from_pressure(max(p,1.0))
    return (1+A[i])*rho*C**2 + max(p,1.0)/(GAMMAS[i]-1)

def integrate_star(pc, rtol=2e-8, max_step=2e4):
    rho_c,_=rho_from_pressure(pc)
    r0=1.0
    epsc=energy_density(pc)
    m0=4*np.pi*r0**3*epsc/(3*C**2)
    den=1-2*G*m0/(r0*C**2)
    dp0=-G*(epsc+pc)*(m0+4*np.pi*r0**3*pc/C**2)/(C**2*r0**2*den)
    p0=pc+0.5*dp0*r0**2
    def rhs(r,y):
        m,p=y; pe=max(float(p),1.0); eps=energy_density(pe)
        fac=max(1-2*G*m/(r*C**2),1e-12)
        return [4*np.pi*r*r*eps/C**2, -G*(eps+pe)*(m+4*np.pi*r**3*pe/C**2)/(C**2*r*r*fac)]
    def surface(r,y): return y[1]
    surface.terminal=True; surface.direction=-1
    sol=solve_ivp(rhs,(r0,3e7),[m0,p0],events=surface,rtol=rtol,atol=[1e15,1e14],max_step=max_step,method='DOP853')
    if not sol.success or len(sol.t_events[0])==0: return None
    R=sol.t_events[0][0]; M=sol.y_events[0][0][0]
    return {'central_pressure_dyn_cm2':float(pc),'central_density_g_cm3':float(rho_c),'mass_solar':float(M/MSUN),'radius_km':float(R/KM)}

def main():
    out=Path(__file__).resolve().parent/'outputs'; out.mkdir(exist_ok=True)
    pcs=np.logspace(34.0,36.5,70)
    rows=[x for pc in pcs if (x:=integrate_star(pc)) is not None]
    if len(rows)<10: raise RuntimeError(f'Only {len(rows)} successful models')
    with (out/'sly_mass_radius_results.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    masses=np.array([r['mass_solar'] for r in rows]); radii=np.array([r['radius_km'] for r in rows])
    imax=int(np.argmax(masses)); i14=int(np.argmin(np.abs(masses-1.4)))
    fig,ax=plt.subplots(figsize=(8,5.5)); ax.plot(radii,masses,'o-',ms=3,lw=1.5,label='SLy piecewise-polytrope approximation')
    ax.scatter([radii[imax]],[masses[imax]],marker='*',s=130,label=f'Max sampled: {masses[imax]:.3f} M☉')
    ax.scatter([radii[i14]],[masses[i14]],marker='s',s=45,label=f'Nearest 1.4 M☉ model: R={radii[i14]:.2f} km')
    ax.set_xlabel('Radius (km)'); ax.set_ylabel('Gravitational mass (solar masses)'); ax.set_title('Neutron-star mass–radius relation — SLy piecewise-polytropic EOS')
    ax.grid(True,alpha=.3); ax.legend(); fig.tight_layout(); fig.savefig(out/'sly_mass_radius_curve.png',dpi=200); plt.close(fig)
    summary=(f'SLy piecewise-polytropic TOV run\nEOS parameters from Read et al. (2009): log10(p1 [dyn/cm^2])=34.384, '
             f'Gamma1=3.005, Gamma2=2.988, Gamma3=2.851; crust/core match rho0=10^14.3, rho1=10^14.7, rho2=10^15 g/cm^3.\n'
             f'Successful models: {len(rows)}/{len(pcs)}\nMaximum mass on sampled grid: {masses[imax]:.5f} Msun at R={radii[imax]:.4f} km.\n'
             f'Nearest sampled 1.4 Msun model: M={masses[i14]:.5f} Msun, R={radii[i14]:.4f} km.\n'
             'Caveat: the outer crust is represented by a matched low-density polytrope, not a tabulated BPS/SLy4 crust. This is a literature-calibrated approximation, not a full microphysical EOS reproduction. Compare the trend to Read et al. (2009), who report approximately Mmax=2.06 Msun and R1.4=11.4 km for SLy.\n')
    (out/'sly_run_summary.txt').write_text(summary)
    print(summary)
if __name__=='__main__': main()
