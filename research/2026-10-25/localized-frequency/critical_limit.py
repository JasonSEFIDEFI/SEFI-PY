"""Formal near-lambda=16 constrained-curve limit: numeric predictions."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_bvp,simpson
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize_scalar,brentq
HERE=Path(__file__).resolve().parent

def solve_limit(beta,seed=None,X=40,tol=1e-9):
    x=np.linspace(0,X,1001)
    if seed is None:
        data=np.load(HERE/'critical_recovery/lambda15.75_beta0.125.npz')
        f=CubicHermiteSpline(data['x'],data['y'][0],data['y'][1])(x)
        s=CubicHermiteSpline(data['x'],data['y'][2],data['y'][3])(x)
        theta=np.arctan2(f/np.sqrt(.75),s/np.sqrt(3/15.75))
        theta[x>32]=0
        theta=np.clip(theta,0,np.pi/2)
        y=np.array([theta,np.gradient(theta,x)])
    else:y=seed.sol(x)
    y[1,0]=0
    def fun(x,y):
        t,v=y;sn=np.sin(t);cs=np.cos(t)
        a=3/16*(1+3*cs*cs);ap=-9/8*sn*cs
        up=3/4*sn*cs*(beta-2*sn*sn)
        return np.array([v,up/(a*beta)-ap*v*v/(2*a)])
    def bc(a,b):return np.array([a[1],b[0]])
    S=np.zeros((2,2));S[1,1]=-2
    sol=solve_bvp(fun,bc,x,y,S=S,tol=tol,bc_tol=1e-12,max_nodes=20000)
    x=np.linspace(0,X,20001);t,v=sol.sol(x);z=x/np.sqrt(beta)
    rho=.75*np.sin(t)**2
    a=3/16*(1+3*np.cos(t)**2)
    K=simpson(4*np.pi*z*z*a*beta*v*v/2,x=z)
    R4=simpson(4*np.pi*z*z*rho*rho,x=z)
    I=simpson(4*np.pi*z*z*rho,x=z)
    Qbar=np.sqrt(1.75)*I
    r=dict(beta=beta,status=int(sol.status),nodes=len(sol.x),X=X,tolerance=tol,
        max_rms_residual=float(max(sol.rms_residuals)),central_theta=float(t[0]),
        scaled_charge=float(Qbar),scaled_binding_difference=float((K-2*R4/3)/Qbar),
        scope='Prediction of a formal constrained-curve limit, not a proved singular-limit theorem')
    assert sol.status==0 and .01<t[0]<np.pi/2,r
    return r,sol

def main():
    out=HERE/'critical_limit_refined.json'
    if out.exists():raise FileExistsError(out)
    result=dict(scope=__doc__,runs=[])
    cache=[]
    def evaluate(beta):
        near=min(cache,key=lambda r:abs(r[0]['beta']-beta))[1] if cache else None
        r,sol=solve_limit(beta,near)
        result['runs'].append(r);cache.append((r,sol))
        out.write_text(json.dumps(result,indent=2)+'\n')
        return r
    for beta in [.125,.15,.175,.2,.225,.25,.275,.3,.1,.09,.08,.075,.065,.06,.055,.05]:evaluate(beta)
    opt=minimize_scalar(lambda b:evaluate(b)['scaled_charge'],bounds=(.1,.15),
        method='bounded',options=dict(xatol=2e-8))
    result['fold_prediction']=evaluate(float(opt.x))
    root=brentq(lambda b:evaluate(b)['scaled_binding_difference'],.25,.275,xtol=1e-10)
    result['binding_prediction']=evaluate(root)
    xx=np.linspace(0,40,1601)
    fold_sol=min(cache,key=lambda r:abs(r[0]['beta']-float(opt.x)))[1]
    binding_sol=min(cache,key=lambda r:abs(r[0]['beta']-root))[1]
    np.savez_compressed(HERE/'critical_limit_profiles.npz',x=xx,
        fold=fold_sol.sol(xx),binding=binding_sol.sol(xx),
        beta_fold=float(opt.x),beta_binding=root)
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='runs'}),flush=True)
if __name__=='__main__':main()

