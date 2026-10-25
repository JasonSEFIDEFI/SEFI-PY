"""Independent continuum collocation check of selected finite-volume profiles."""
from pathlib import Path
import json, platform
import numpy as np
import scipy
from scipy.integrate import solve_bvp, cumulative_trapezoid, simpson
HERE=Path(__file__).resolve().parent

def validate(folder,Q,tol=1e-7,seed_Q=None):
    report=json.loads((HERE/folder/'results.json').read_text())
    seed_Q=Q if seed_Q is None else seed_Q
    meta=next(r for r in report['runs'] if r['Q']==seed_Q)
    data=np.load(HERE/folder/('Q'+format(seed_Q,'.10g')+'.npz'))
    L=meta['L'];lam=meta['lambda_sigma'];vac=np.sqrt(3/lam)
    x=np.linspace(0,L,round(L/.2)+1)
    f=np.interp(x,data['r'],data['f'],right=0)
    s=np.interp(x,data['r'],data['s'],right=vac)
    J=cumulative_trapezoid(4*np.pi*x*x*f*f,x,initial=0)
    y=np.array([f,np.gradient(f,x),s,np.gradient(s,x),J])
    y[1,0]=0;y[3,0]=0
    def fun(r,y,p):
        f,fp,s,sp,J=y;w=p[0]
        return np.array([fp,(1+f*f+4*s*s-w*w)*f,
                         sp,(-3+4*f*f+lam*s*s)*s,4*np.pi*r*r*f*f])
    def bc(a,b,p):
        return np.array([a[1],a[3],a[4],b[0],b[2]-vac,p[0]*b[4]-Q])
    S=np.zeros((5,5));S[1,1]=-2;S[3,3]=-2
    sol=solve_bvp(fun,bc,x,y,p=[meta['omega']],S=S,tol=tol,bc_tol=1e-8,max_nodes=30000)
    xx=np.linspace(0,L,20001);f,fp,s,sp,J=sol.sol(xx);w=sol.p[0]
    rho=f*f;eta=s*s
    potential=rho*rho/4+rho/2+2*rho*eta+lam*(eta-3/lam)**2/4
    I=simpson(4*np.pi*xx*xx*f*f,x=xx)
    K=simpson(4*np.pi*xx*xx*(fp*fp+sp*sp)/2,x=xx)
    U=simpson(4*np.pi*xx*xx*potential,x=xx)
    C=w*w*I/2;E=K+U+C
    r=dict(source_folder=folder,seed_Q=seed_Q,Q=Q,L=L,lambda_sigma=lam,
        status=int(sol.status),message=sol.message,tolerance=tol,nodes=len(sol.x),
        max_rms_residual=float(max(sol.rms_residuals)),
        boundary_residual=float(max(abs(bc(sol.y[:,0],sol.y[:,-1],sol.p)))),
        charge_relative_error=float(abs(w*I-Q)/Q),
        energy=float(E),energy_per_charge=float(E/Q),omega=float(w),
        finite_volume_energy_per_charge=meta['energy_per_charge'] if Q==seed_Q else None,
        relative_energy_difference=float((E-meta['energy'])/E) if Q==seed_Q else None,
        virial_relative=float((K+3*U-3*C)/E))
    return r

if __name__=='__main__':
    out=HERE/'independent_bvp_results.json'
    if out.exists():raise FileExistsError(out)
    result=dict(scope=__doc__,python=platform.python_version(),scipy=scipy.__version__,runs=[])
    for folder,Q in [('h0p1_L80',100000.),('h0p1_L80',5000.),
                     ('h0p1_L80',1250.),('threshold_h0p2_L100',1200.),
                     ('threshold_h0p2_L100',1100.)]:
        r=validate(folder,Q);result['runs'].append(r)
        out.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(r),flush=True)

