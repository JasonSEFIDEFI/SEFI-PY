"""Explicit compact spherical charge trial, not a vortex or solved soliton.

Uses the existing physical-potential family over its zero-charge vacuum.
Quintic wall makes the radial integrands polynomial of degree <=22;
12-point Gauss-Legendre integrates them exactly in exact arithmetic.
"""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent

def trial(Q,lam,R,width,amplitude,order=12):
    if R<=width or width<=0:raise ValueError('Require R > width > 0')
    nodes,weights=np.polynomial.legendre.leggauss(order);x=(nodes+1)/2
    r=R-width+2*width*x;weights=weights*width
    t=1-10*x**3+15*x**4-6*x**5
    dt=-30*x*x*(1-x)**2/(2*width)
    f=amplitude*t;s=np.sqrt(3/lam)*(1-t)
    rho=f*f;eta=s*s;vol=4*np.pi*r*r*weights
    inner_volume=4*np.pi*(R-width)**3/3
    I=float(np.sum(vol*rho)+inner_volume*amplitude**2)
    K=float(.5*np.sum(vol*(amplitude**2+3/lam)*dt**2))
    def potential(rho,eta):return rho*rho/4+rho/2+2*rho*eta+lam*(eta-3/lam)**2/4
    U=float(np.sum(vol*potential(rho,eta))+inner_volume*potential(amplitude**2,0))
    energy=K+U+Q*Q/(2*I)
    return dict(Q=Q,lam=lam,R=R,width=width,amplitude=amplitude,I=I,omega=Q/I,
                gradient_energy=K,potential_energy=U,energy=energy,energy_per_charge=energy/Q,
                free_phi_mass=float(np.sqrt(1+12/lam)))

if __name__=='__main__':
    output=[]
    for lam in [25.,12.]:
        for Q in [1e2,1e3,1e4,1e5,1e6,1e7]:
            best=None
            for a in [.75,.85,.93,1.05,1.2]:
                center=(3*Q/(4*np.pi*np.sqrt(1+12/lam)*a*a))**(1/3)
                for scale in np.linspace(.8,1.4,19):
                    R=center*scale
                    for width in [.3,.5,.75,1.,1.5,2.]:
                        if R<=width:continue
                        result=trial(Q,lam,R,width,a)
                        if best is None or result['energy']<best['energy']:best=result
            refined=trial(Q,lam,best['R'],best['width'],best['amplitude'],24)
            best['quadrature_refinement_relative_difference']=abs(refined['energy']-best['energy'])/abs(refined['energy'])
            best['below_free_threshold']=best['energy_per_charge']<best['free_phi_mass']
            if lam==25:assert best['energy_per_charge']>=best['free_phi_mass']-1e-10
            output.append(best);print(json.dumps(best),flush=True)
    (ROOT/'vacuum_trial.json').write_text(json.dumps(dict(scope='Compact trial configurations; no stationary PDE solution or stability/physical-matter identification',cases=output),indent=2)+'\n')
