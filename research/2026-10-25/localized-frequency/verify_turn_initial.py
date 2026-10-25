"""Tighter fold, equal-charge branch pair, and coupling binding crossings."""
from pathlib import Path
from types import SimpleNamespace
import json
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize_scalar,brentq
from frequency_scan import solve,save_profile
HERE=Path(__file__).resolve().parent

def seed_from(record,folder):
    name='lambda'+format(record['lambda_sigma'],'.6g')+'_beta'+format(record['beta'],'.6g')+'.npz'
    data=np.load(HERE/folder/name)
    f=CubicHermiteSpline(data['x'],data['y'][0],data['y'][1])
    s=CubicHermiteSpline(data['x'],data['y'][2],data['y'][3])
    def evaluate(x):return np.array([f(x),f(x,1),s(x),s(x,1)])
    return record,SimpleNamespace(sol=evaluate)

def main():
    out=HERE/'verification';out.mkdir(exist_ok=False)
    scan=json.loads((HERE/'scan_X32/results.json').read_text())
    result=dict(scope=__doc__,evaluations=[],binding_crossings=[],paired_states=[])
    def save():(out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    def evaluate(lam,beta,X=40,tol=1e-9):
        candidates=[r for r in scan['runs'] if r['accepted'] and r['lambda_sigma']==lam]
        near=min(candidates,key=lambda r:abs(r['beta']-beta))
        r,sol=solve(lam,beta,seed_from(near,'scan_X32'),X=X,tol=tol)
        result['evaluations'].append(r);save()
        if not r['accepted']:raise RuntimeError(str(r))
        return r,sol
    opt=minimize_scalar(lambda beta:evaluate(12.,beta)[0]['Q'],
        bounds=(.1,.15),method='bounded',options=dict(xatol=2e-7))
    fold,sol=evaluate(12.,float(opt.x));result['refined_fold']=fold
    save_profile(out/'refined_fold_lambda12.npz',fold,sol);save()
    for label,interval in [('high_frequency',(.025,.12)),('low_frequency',(.13,.275))]:
        beta=brentq(lambda b:evaluate(12.,b)[0]['Q']-1000,*interval,xtol=1e-10)
        r,sol=evaluate(12.,beta);r['branch_label']=label
        result['paired_states'].append(r)
        save_profile(out/(label+'_Q1000.npz'),r,sol);save()
    for lam in sorted(set(r['lambda_sigma'] for r in scan['runs'])):
        beta=brentq(lambda b:evaluate(lam,b,X=40,tol=1e-8)[0]['binding_difference'],
            .2,.35,xtol=1e-8)
        r,sol=evaluate(lam,beta,X=40,tol=1e-8)
        result['binding_crossings'].append(r)
        save_profile(out/('binding_lambda'+format(lam,'.6g')+'.npz'),r,sol)
        save();print(json.dumps(dict(lambda_sigma=lam,beta=beta,Q_binding=r['Q'])),flush=True)
    save()
    print(json.dumps(dict(refined_fold=result['refined_fold'],
                         paired_states=result['paired_states'])),flush=True)
if __name__=='__main__':main()

