"""Check reduced predictions with a larger domain and tighter tolerance."""
from pathlib import Path
from types import SimpleNamespace
import json
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from critical_limit import solve_limit
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    out=HERE/'critical_limit_validation.json'
    if out.exists():raise FileExistsError(out)
    meta=json.loads((HERE/'critical_limit_refined.json').read_text())
    data=np.load(HERE/'critical_limit_profiles.npz')
    result=dict(scope=__doc__,runs=[])
    for key in ['fold','binding']:
        y=data[key];curve=CubicHermiteSpline(data['x'],y[0],y[1])
        def evaluate(x):
            clipped=np.minimum(x,data['x'][-1])
            t=curve(clipped);v=curve(clipped,1)
            t[x>data['x'][-1]]=0;v[x>data['x'][-1]]=0
            return np.array([t,v])
        baseline=meta[key+'_prediction']
        r,sol=solve_limit(baseline['beta'],SimpleNamespace(sol=evaluate),X=48,tol=1e-10)
        r['charge_difference_from_X40']=r['scaled_charge']-baseline['scaled_charge']
        r['scaled_binding_difference_from_X40']=r['scaled_binding_difference']-baseline['scaled_binding_difference']
        result['runs'].append(r)
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

