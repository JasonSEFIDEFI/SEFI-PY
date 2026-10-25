"""Compare fixed-Q radial amplitude curvatures of the two Q=1000 states."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.sparse.linalg import eigsh
HERE=Path(__file__).resolve().parent
parent=HERE/'reference'
sys.path.insert(0,str(parent))
from vacuum_spherical import Radial
from continue_branch import polish
from spherical_stability import operators,check_hessian

def main():
    out=HERE/'paired_hessian_refined.json'
    if out.exists():raise FileExistsError(out)
    report=dict(scope=__doc__,hessian_directional_relative_error=check_hessian(),runs=[])
    for label in ['high_frequency','low_frequency']:
        data=np.load(HERE/'verification_refined'/(label+'_Q1000.npz'))
        ffun=CubicHermiteSpline(data['x'],data['y'][0],data['y'][1])
        sfun=CubicHermiteSpline(data['x'],data['y'][2],data['y'][3])
        kappa=float(data['kappa'])
        for L,h in [(120.,.4),(120.,.2),(120.,.1),(160.,.2)]:
            model=Radial(Q=1000,lam=12,L=L,h=h)
            f=ffun(model.r*kappa);s=sfun(model.r*kappa)
            f,s,it=polish(model,f,s,1e-9)
            gf,gs,w=model.gradient(f,s)
            res=float(max(abs(gf).max(),abs(gs).max()))
            assert res<1e-9,res
            op,_,_=operators(model,f,s,0)
            vals,vecs=eigsh(op,k=3,which='SA',tol=1e-9,maxiter=30000,
                v0=np.random.default_rng(12).normal(size=op.shape[0]))
            r=dict(branch_label=label,L=L,h=h,Q=1000,omega=w,
                energy=model.energy(f,s),stationary_residual=res,
                eigenvalues=vals.tolist(),
                eigenpair_residuals=[float(np.linalg.norm(op@vecs[:,i]-vals[i]*vecs[:,i])) for i in range(3)])
            report['runs'].append(r)
            out.write_text(json.dumps(report,indent=2)+'\n')
            print(json.dumps(r),flush=True)
if __name__=='__main__':main()

