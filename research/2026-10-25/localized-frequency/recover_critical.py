"""Small frequency steps near lambda=16; failed initial attempts are retained."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize_scalar
from frequency_scan import solve,save_profile
from verify_turn import seed_from
HERE=Path(__file__).resolve().parent

def main():
    out=HERE/'critical_recovery';out.mkdir(exist_ok=False)
    initial=json.loads((HERE/'scan_X32/results.json').read_text())
    report=dict(scope=__doc__,runs=[],fold_evaluations=[],folds=[])
    def save():(out/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    targets=[.25,.225,.2,.175,.15,.125,.1,.075,.05,.025,.0125,.00625]
    for lam in [15.5,15.75,15.9]:
        anchor=next(r for r in initial['runs'] if r['lambda_sigma']==lam and r['beta']==.275 and r['accepted'])
        seed=seed_from(anchor,'scan_X32');cache=[seed]
        for target in targets:
            n=max(1,int(np.ceil((seed[0]['beta']-target)/(.00125 if target<.05 else .005))))
            for b in np.linspace(seed[0]['beta'],target,n+1)[1:]:
                r,sol=solve(lam,float(b),seed,X=32,tol=1e-7)
                report['runs'].append(r);save()
                if not r['accepted']:raise RuntimeError(str(r))
                seed=(r,sol);cache.append(seed)
            save_profile(out/('lambda'+format(lam,'.6g')+'_beta'+format(target,'.6g')+'.npz'),r,sol)
            print(json.dumps(dict(lambda_sigma=lam,beta=target,Q=r['Q'])),flush=True)
        def objective(beta):
            near=min(cache,key=lambda z:abs(z[0]['beta']-beta))
            r,sol=solve(lam,beta,near,X=40,tol=1e-8)
            report['fold_evaluations'].append(r);save()
            if not r['accepted']:raise RuntimeError(str(r))
            cache.append((r,sol));return r['Q']
        opt=minimize_scalar(objective,bounds=(.11,.14),method='bounded',options=dict(xatol=3e-7))
        best=min(cache,key=lambda z:z[0]['Q'])
        report['folds'].append(dict(lambda_sigma=lam,beta=float(opt.x),Q=float(opt.fun),profile=best[0]))
        save_profile(out/('fold_lambda'+format(lam,'.6g')+'.npz'),best[0],best[1])
        save()
    save()
if __name__=='__main__':main()

