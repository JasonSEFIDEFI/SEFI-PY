"""Radial fixed-charge energy-curvature crossing near full-model charge folds."""
import json,time
import numpy as np
from scipy.linalg import eigh
from session import ROOT,reduced_seed,full_center,save

def matrix(sol,delta,alpha,L,h):
    n=round(L/h);h=L/n;r=np.arange(1,n)*h
    f,fp,s,sp=sol.sol(r);lam=16-delta
    omega2=1+12/lam-delta*alpha
    off=np.full(len(r)-1,-1/h**2)
    K=np.diag(np.full(len(r),2/h**2))+np.diag(off,1)+np.diag(off,-1)
    Vf=(1+3*f*f+4*s*s-omega2)/delta
    Vs=(-3+4*f*f+3*lam*s*s)/delta
    C=8*f*s/delta
    H=np.block([[K+np.diag(Vf),np.diag(C)],
                [np.diag(C),K+np.diag(Vs)]])
    J=h*np.sum(r*r*f*f);u=np.r_[np.sqrt(h)*r*f,np.zeros(len(r))]
    H+=4*omega2/(delta*J)*np.outer(u,u)
    return H,r

def main():
    import session
    session.SOLVE_FILE='radial_profile_solves.json'
    folds=json.loads((ROOT/'folds.json').read_text())
    seed=reduced_seed();rows=[]
    for delta in [.25,.03125]:
        fold=next(r for r in folds if r['delta']==delta)
        c0=fold['center'];profiles={}
        for offset in [-.025,0.,.025]:
            sol,row=full_center(delta,c0+offset,seed)
            profiles[offset]=(sol,row)
        # Finite derivative of the stationary family for the fold's null direction.
        eps=1e-4
        minus,mrow=full_center(delta,c0-eps,seed)
        plus,prow=full_center(delta,c0+eps,seed)
        dq=(prow['q_scaled']-mrow['q_scaled'])/(2*eps)
        for offset,(sol,meta) in profiles.items():
            for h in [.5,.25]:
                start=time.perf_counter();H,r=matrix(sol,delta,meta['alpha'],240,h)
                vals,vec=eigh(H,subset_by_index=[0,1],driver='evr',check_finite=False)
                row=dict(delta=delta,center_offset=offset,center=meta['center'],
                    h=h,L=240,scaled_lowest_eigenvalues=vals.tolist(),
                    eigen_residual=float(max(np.linalg.norm(H@vec[:,i]-vals[i]*vec[:,i]) for i in range(2))),
                    seconds=time.perf_counter()-start)
                if offset==0:
                    tangent=(plus.sol(r)-minus.sol(r))/(2*eps)
                    u=np.r_[np.sqrt(h)*r*tangent[0],np.sqrt(h)*r*tangent[2]]
                    rhs=np.r_[np.sqrt(h)*r*sol.sol(r)[0],np.zeros(len(r))]
                    rhs*=2*meta['omega']*dq/(delta*meta['I_scaled'])
                    row.update(charge_derivative_c=dq,
                        tangent_rayleigh=float(u@(H@u)/(u@u)),
                        family_identity_residual_per_tangent=float(np.linalg.norm(H@u-rhs)/np.linalg.norm(u)),
                        lowest_vector_tangent_overlap=float(abs(vec[:,0]@u)/np.linalg.norm(u)))
                rows.append(row);save('radial_results.json',rows)
                print(json.dumps(row),flush=True)
    # Keep this process's solves separate from the independent minimum-charge sweep.
    from session import SOLVES
    save('radial_profile_solves.json',SOLVES)

if __name__=='__main__':main()
