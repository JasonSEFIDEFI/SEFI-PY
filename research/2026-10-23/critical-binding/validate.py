"""Domain/tolerance checks and independent full two-field comparison."""
import json,time,platform
import numpy as np
from scipy.integrate import solve_bvp,simpson
import scipy
from critical import solve,ROOT

checks=[]
warm=None
for c in [.5,.7,.9,1.1,1.3,1.45]:
    warm,row=solve(c,guess=warm)
    assert row['admissible'],row
seed=warm
for center in [1.1,1.45]:
    ref=seed
    for L,n,tol in [(160,1000,1e-6),(240,1600,1e-8),(320,2400,1e-9)]:
        sol,row=solve(center,L,n,tol,ref);checks.append(row);ref=sol
        print('reduced',json.dumps(row),flush=True)

center=1.45
red,reference=solve(center,240,1600,1e-9,seed)
assert reference['admissible'],reference
full=[]
for delta in [1.,.5,.25,.125,.0625,.03125,.015625,.0078125]:
    lam=16-delta;vac=np.sqrt(3/lam);m=np.sqrt(1+12/lam)
    x=np.linspace(0,240,1600);t,v=red.sol(x)
    y=np.array([np.sqrt(3)/2*np.sin(t),np.sqrt(3)/2*np.cos(t)*v,
                vac*np.cos(t),-vac*np.sin(t)*v])
    def fun(x,y,p):
        f,fp,s,ss=y;a=p[0]
        return np.array([fp,((1-m*m+delta*a+f*f+4*s*s)*f)/delta,
                         ss,((-3+4*f*f+lam*s*s)*s)/delta])
    def bc(ya,yb,p):
        return np.array([ya[0]-np.sqrt(3)/2*np.sin(center),ya[1],ya[3],yb[0],yb[2]-vac])
    S=np.diag([0.,-2.,0.,-2.]);start=time.perf_counter()
    sol=solve_bvp(fun,bc,x,y,p=red.p,S=S,tol=1e-7,max_nodes=40000)
    x=np.linspace(0,240,20001);f,fp,s,ss=sol.sol(x);a=float(sol.p[0])
    omega=np.sqrt(m*m-delta*a)
    assert sol.success and 0<a<(m*m-(1+3/np.sqrt(lam)))/delta
    assert f.min()>-1e-7 and s.min()>-1e-7 and np.max(fp)<1e-7
    I=4*np.pi*simpson(x*x*f*f,x=x);q=omega*I
    # Compute subtracted potential through the exact square identity to avoid cancellation.
    diff=lam/4*(s*s+(4*f*f-3)/lam)**2+(lam-16)*f**4/(4*lam)
    action=4*np.pi*simpson(x*x*(.5*(fp*fp+ss*ss)+diff/delta+a*f*f/2),x=x)
    binding=action-a*q/(m+omega)
    row=dict(delta=delta,lam=lam,success=bool(sol.success),alpha=a,q_scaled=q,
             binding_scaled=binding,max_rms_residual=float(max(sol.rms_residuals)),
             relative_q_error=q/reference['q_scaled']-1,
             seconds=time.perf_counter()-start,nodes=len(sol.x),message=sol.message)
    full.append(row);print('full',json.dumps(row),flush=True)
    np.savez_compressed(ROOT/f'full_delta_{delta}.npz',rho=x,f=f,s=s,alpha=a)
    (ROOT/'validation.json').write_text(json.dumps(dict(reduced=checks,reference=reference,full=full,
         environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)),indent=2))
