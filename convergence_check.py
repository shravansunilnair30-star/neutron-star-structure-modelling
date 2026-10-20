
import sys, csv
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import neutron_star_tov as tov
outdir = Path(__file__).resolve().parent / 'outputs'
outdir.mkdir(exist_ok=True)

def solve(p_c, rtol=1e-7, max_step=1e4):
    k, gamma = tov.K, tov.GAMMA
    rho_c = (p_c/k)**(1/gamma)
    r0=1.0
    eps_c = rho_c*tov.C**2 + p_c/(gamma-1)
    m0=(4/3)*np.pi*r0**3*eps_c/tov.C**2
    denom=1-2*tov.G*m0/(r0*tov.C**2)
    dpdr=-tov.G*(eps_c+p_c)*(m0+4*np.pi*r0**3*p_c/tov.C**2)/(tov.C**2*r0**2*denom)
    p0=p_c+0.5*dpdr*r0**2
    def rhs(r,y):
        m,p=y; pe=max(p,1.0); rho=(pe/k)**(1/gamma); eps=rho*tov.C**2+pe/(gamma-1)
        fac=max(1-2*tov.G*m/(r*tov.C**2),1e-12)
        return [4*np.pi*r*r*eps/tov.C**2,
                -tov.G*(eps+pe)*(m+4*np.pi*r**3*pe/tov.C**2)/(tov.C**2*r*r*fac)]
    def ev(r,y): return y[1]
    ev.terminal=True; ev.direction=-1
    sol=solve_ivp(rhs,(r0,5e7),[m0,p0],events=ev,rtol=rtol,atol=[1e15,1e15],max_step=max_step)
    if not sol.success or not len(sol.t_events[0]): return None
    return sol.y_events[0][0][0]/tov.M_SUN, sol.t_events[0][0]/tov.KM

pcs=np.logspace(33,36,40)
configs=[("baseline",1e-7,1e4),("tighter_tol",1e-9,1e4),("smaller_step",1e-7,5e3),("both_tighter",1e-9,5e3)]
rows=[]
for name,rtol,step in configs:
    vals=[]
    for pc in pcs:
        v=solve(pc,rtol,step)
        if v: vals.append((pc,*v))
    if vals:
        arr=np.array(vals)
        i=np.argmax(arr[:,1])
        rows.append({"configuration":name,"rtol":rtol,"max_step_cm":step,"successful_models":len(vals),
                     "max_sampled_mass_Msun":arr[i,1],"radius_at_max_km":arr[i,2],
                     "central_pressure_at_max_dyn_cm2":arr[i,0]})
with open(outdir/"convergence_test.csv","w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)
with open(outdir/"convergence_test_summary.txt","w",encoding="utf-8") as f:
    f.write("Convergence check for the simplified polytropic TOV model\\n\\n")
    f.write("Configurations vary ODE relative tolerance and maximum radial step.\\n")
    f.write("This checks numerical sensitivity, not physical accuracy of the EOS.\\n\\n")
    for r in rows: f.write(str(r)+"\\n")
    if len(rows)>1:
        base=rows[0]
        for r in rows[1:]:
            dm=abs(r["max_sampled_mass_Msun"]-base["max_sampled_mass_Msun"])/abs(base["max_sampled_mass_Msun"])*100
            dr=abs(r["radius_at_max_km"]-base["radius_at_max_km"])/abs(base["radius_at_max_km"])*100
            f.write(f"Relative change vs baseline for {r['configuration']}: mass {dm:.6g}%, radius {dr:.6g}%\\n")
print(rows)
