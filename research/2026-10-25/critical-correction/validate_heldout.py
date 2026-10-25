"""Larger-domain/tighter-tolerance check of final held-out stationary profiles."""
from pathlib import Path
from types import SimpleNamespace
import json
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from full_solver import solve,save_profile
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    out=HERE/'heldout_validation';out.mkdir(exist_ok=False)
    original=json.loads((HERE/'heldout/results.json').read_text())
    result=dict(scope=__doc__,runs=[])
    for test in original['tests']:
        tag=test['kind']+'_lambda'+format(test['lambda_sigma'],'.6g')
        data=np.load(HERE/'heldout'/(tag+'.npz'))
        f=CubicHermiteSpline(data['x'],data['y'][0],data['y'][1])
        s=CubicHermiteSpline(data['x'],data['y'][2],data['y'][3])
        def evaluate(x):return np.array([f(x),f(x,1),s(x),s(x,1)])
        seed=(test['full_model_record'],SimpleNamespace(sol=evaluate))
        r,sol=solve(test['lambda_sigma'],test['beta'],seed,X=48,tol=1e-10)
        if not r['accepted']:raise RuntimeError(str(r))
        r['kind']=test['kind']
        r['scaled_charge_change']=r['scaled_charge']-test['scaled_actual']
        r['relative_scaled_charge_change']=r['scaled_charge_change']/r['scaled_charge']
        r['check_scope']='Same beta as held-out optimum/root; tests solver/domain changes, not a reoptimization.'
        result['runs'].append(r)
        save_profile(out/(tag+'.npz'),r,sol)
        (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(r),flush=True)

