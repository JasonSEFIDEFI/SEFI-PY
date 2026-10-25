"""Threshold, first-law and selected energy-Hessian checks; exploratory."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq
from scipy.sparse.linalg import eigsh
from independent_bvp import validate
from vacuum_spherical import Radial
from spherical_stability import operators,check_hessian
HERE=Path(__file__).resolve().parent

if __name__=='__main__':
    out=HERE/'additional_checks.json'
    if out.exists():raise FileExistsError(out)
    result=dict(scope=__doc__,hessian_directional_relative_error=check_hessian(),
                threshold_runs=[],first_law=[],spectra=[])
    def save():out.write_text(json.dumps(result,indent=2)+'\n')
    folder='threshold_h0p05_L120'
    def objective(Q,tol=1e-8):
        r=validate(folder,Q,tol=tol,seed_Q=1225.)
        assert r['status']==0,r
        result['threshold_runs'].append(r);save()
        return r['energy_per_charge']-np.sqrt(2)
    qstar=brentq(objective,1225.,1230.,xtol=1e-5)
    refined=validate(folder,qstar,tol=1e-9,seed_Q=1225.)
    result['binding_crossing']=dict(Q=qstar,refined_check=refined,
        scope='Numerical E(Q)/Q=sqrt(2) crossing on selected spherical branch, finite L=120; not minimum existence charge')
    save()
    for Q in [1250.,5000.,100000.]:
        folder0='h0p1_L80'
        vals=[]
        for fraction in [.001,.0005]:
            dQ=Q*fraction
            lo=validate(folder0,Q-dQ,tol=1e-8,seed_Q=Q)
            hi=validate(folder0,Q+dQ,tol=1e-8,seed_Q=Q)
            mid=validate(folder0,Q,tol=1e-8)
            assert all(r['status']==0 for r in [lo,hi,mid])
            deriv=(hi['energy']-lo['energy'])/(2*dQ)
            vals.append(dict(relative_step=fraction,dEdQ=deriv,omega=mid['omega'],
                absolute_error=abs(deriv-mid['omega']),runs=[lo,hi,mid]))
        result['first_law'].append(dict(Q=Q,checks=vals));save()
    for folder0,Q in [('h0p1_L80',1250.),('threshold_h0p1_L100',1200.),
                      ('threshold_h0p1_L100',1000.)]:
        report=json.loads((HERE/folder0/'results.json').read_text())
        meta=next(r for r in report['runs'] if r['Q']==Q)
        data=np.load(HERE/folder0/('Q'+format(Q,'.10g')+'.npz'))
        model=Radial(Q=Q,lam=12.,L=meta['L'],h=meta['h'])
        entry=dict(folder=folder0,Q=Q,sectors=[])
        for ell in [0,1,2]:
            sector=dict(ell=ell)
            for name,op in zip(['amplitude','phi_phase','sigma_phase'],
                               operators(model,data['f'],data['s'],ell)):
                vals,vecs=eigsh(op,k=2,which='SA',tol=1e-8,maxiter=30000,
                    v0=np.random.default_rng(12).normal(size=op.shape[0]))
                sector[name]=dict(eigenvalues=vals.tolist(),
                    residuals=[float(np.linalg.norm(op@vecs[:,i]-vals[i]*vecs[:,i])) for i in range(2)])
            entry['sectors'].append(sector)
        result['spectra'].append(entry);save()
        print(json.dumps(entry),flush=True)
    save()
    print(json.dumps(result['binding_crossing']),flush=True)

