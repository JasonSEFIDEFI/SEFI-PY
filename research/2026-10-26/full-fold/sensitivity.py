"""Linear BVP family tangent and matching-domain charge-fold certificate."""
import numpy as np
from scipy.integrate import solve_bvp,simpson
from session import ROOT,reduced_seed,full_center,save
import session,json

def tangent(sol,meta):
    d=meta['delta'];lam=16-d;a=meta['alpha'];w=meta['omega'];L=meta['L']
    x=np.linspace(0,L,2000)
    def fun(r,y,p):
        f,fp,s,sp=sol.sol(r);u,up,v,vp=y
        A=(1-w*w+3*f*f+4*s*s)/d
        D=(-3+4*f*f+3*lam*s*s)/d;C=8*f*s/d
        return np.array([up,A*u+C*v+p[0]*f,vp,C*u+D*v])
    def bc(ya,yb,p):
        return np.array([ya[0]-np.sqrt(3)/2*np.cos(meta['center']),ya[1],ya[3],yb[0],yb[2]])
    ans=solve_bvp(fun,bc,x,np.zeros((4,len(x))),p=[.006],
        S=np.diag([0.,-2.,0.,-2.]),tol=2e-10,max_nodes=40000)
    assert ans.success,ans.message
    r=np.linspace(0,L,24001);f=sol.sol(r)[0];u=ans.sol(r)[0]
    Ic=8*np.pi*simpson(r*r*f*u,x=r)
    wc=-d*float(ans.p[0])/(2*w)
    qc=wc*meta['I_scaled']+w*Ic
    return ans,dict(alpha_c=float(ans.p[0]),I_c=float(Ic),omega_c=wc,q_c=float(qc),
                    bvp_residual=float(max(ans.rms_residuals)),nodes=len(ans.x))

def certificate(sol,meta,tan,deriv,h):
    d=meta['delta'];lam=16-d;L=meta['L'];w=meta['omega']
    n=round(L/h);r=np.arange(1,n)*h
    f,fp,s,sp=sol.sol(r);u,up,v,vp=tan.sol(r)
    y=np.r_[np.sqrt(h)*r*u,np.sqrt(h)*r*v]
    N=len(r);a=y[:N];b=y[N:]
    def kin(z):return (2*z-np.r_[0,z[:-1]]-np.r_[z[1:],0])/h**2
    A=(1-w*w+3*f*f+4*s*s)/d
    D=(-3+4*f*f+3*lam*s*s)/d;C=8*f*s/d
    out=np.r_[kin(a)+A*a+C*b,kin(b)+C*a+D*b]
    z=np.r_[np.sqrt(h)*r*f,np.zeros(N)]
    J=h*np.sum(r*r*f*f)
    out+=4*w*w/(d*J)*z*np.dot(z,y)
    rhs=2*w*deriv['q_c']/(d*meta['I_scaled'])*z
    return dict(h=h,L=L,identity_residual_per_tangent=float(np.linalg.norm(out-rhs)/np.linalg.norm(y)),
                rayleigh=float(y@out/(y@y)),q_c=deriv['q_c'])

def main():
    session.SOLVE_FILE='sensitivity_profile_solves.json'
    data=json.loads((ROOT/'results.json').read_text());seed=reduced_seed();rows=[]
    for d in [.25,.03125]:
        f=next(r for r in data['folds'] if r['delta']==d)
        sol,meta=full_center(d,f['center'],seed,tol=2e-11)
        tan,deriv=tangent(sol,meta)
        # One derivative-based Newton correction to the charge-minimizing center.
        corrected_c=f['center']-deriv['q_c']/f['charge_second_derivative_c']
        sol,meta=full_center(d,corrected_c,seed,tol=2e-11)
        tan,corrected=tangent(sol,meta)
        checks=[certificate(sol,meta,tan,corrected,h) for h in [.5,.25,.125]]
        row=dict(delta=d,original_center=f['center'],corrected_center=corrected_c,
            original_tangent=deriv,corrected_tangent=corrected,checks=checks)
        rows.append(row);save('sensitivity_results.json',rows)
        np.savez_compressed(ROOT/f'tangent_delta_{d}.npz',rho=tan.x,fields_c=tan.y)
        print(json.dumps(row),flush=True)

if __name__=='__main__':main()
