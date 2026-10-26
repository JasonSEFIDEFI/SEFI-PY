"""Static fixed-charge energy curvature for REAL Sigma, complex Phi.

Continuum profiles interpolated to a second-order radial u=r*perturbation grid.
Finite-difference eigenvalues are convergence evidence, not a nonlinear theorem.
"""
import json,time
import numpy as np
from scipy.linalg import eigh,eigh_tridiagonal
from scipy.integrate import simpson
from continue_math import ROOT,central_solve,reduced_alpha,full_alpha

def curvature(sol,delta,a,L,h):
    n=round(L/h);h=L/n;r=np.arange(1,n)*h
    f,fp,s,sp=sol.sol(r);lam=16-delta
    omega2=1+12/lam-delta*a
    J=simpson(np.r_[0.,r*r*f*f,0.],x=np.r_[0.,r,L])
    # J quadrature only sets the global rank-one normalization.
    # Its vector uses the same uniform discrete inner product as the matrix.
    J_discrete=h*np.sum(r*r*f*f)
    assert abs(J/J_discrete-1)<1e-4
    diagonal=np.full(len(r),2/h**2);off=np.full(len(r)-1,-1/h**2)
    K=np.diag(diagonal)+np.diag(off,1)+np.diag(off,-1)
    Vf=(1+3*f*f+4*s*s-omega2)/delta
    Vs=(-3+4*f*f+3*lam*s*s)/delta
    cross=8*f*s/delta
    phaseV=(1+f*f+4*s*s-omega2)/delta
    w=np.r_[np.sqrt(h)*r*f,np.zeros(len(r))]
    records=[]
    for l in [0,1,2]:
        start=time.perf_counter();angular=l*(l+1)/r**2
        H=np.block([[K+np.diag(Vf+angular),np.diag(cross)],
                    [np.diag(cross),K+np.diag(Vs+angular)]])
        if l==0:H+=(4*omega2/(delta*J_discrete))*np.outer(w,w)
        vals,vecs=eigh(H,subset_by_index=[0,2],check_finite=False,driver='evr')
        residual=max(np.linalg.norm(H@vecs[:,i]-vals[i]*vecs[:,i]) for i in range(3))
        pv=eigh_tridiagonal(diagonal+phaseV+angular,off,select='i',select_range=(0,2),check_finite=False)[0]
        records.append(dict(l=l,real_eigenvalues=vals.tolist(),phi_phase_eigenvalues=pv.tolist(),
                            max_eigen_residual=float(residual),seconds=time.perf_counter()-start))
    return dict(h=h,L=L,delta=delta,alpha=a,sectors=records,
                normalization_relative_difference=float(J/J_discrete-1))

def main():
    start=time.perf_counter();previous=None;seeds={}
    for c in [.5,.7,.9,1.1,1.3,1.45]:
        sol,row=central_solve(c,320,1800,1e-9,previous)
        assert row['admissible'],row
        previous=sol
        if c in [.9,1.45]:seeds[c]=sol
    result=dict(scope=__doc__,branches=[])
    for c,seed in seeds.items():
        a=float(seed.p[0]);red,_=reduced_alpha(a,seed)
        sol,row=full_alpha(.03125,a,red)
        branch=dict(reduced_center=c,full_solution=row,grids=[])
        for h in [.5,.25]:
            record=curvature(sol,.03125,a,240,h)
            branch['grids'].append(record)
            print(json.dumps(dict(center=c,**record)),flush=True)
        result['branches'].append(branch)
        (ROOT/'stability_results.json').write_text(json.dumps(result,indent=2))
    result['seconds']=time.perf_counter()-start
    (ROOT/'stability_results.json').write_text(json.dumps(result,indent=2))

if __name__=='__main__':main()
