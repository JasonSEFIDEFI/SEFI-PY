"""Test frozen first-order predictions against closer full-model couplings."""
from pathlib import Path
from types import SimpleNamespace
import json,hashlib
from datetime import datetime,timezone
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize_scalar,brentq
from full_solver import solve,save_profile
HERE=Path(__file__).resolve().parent

def load_seed(kind):
    data=np.load(HERE/('seed_'+kind+'_lambda15p9.npz'))
    f=CubicHermiteSpline(data['x'],data['y'][0],data['y'][1])
    s=CubicHermiteSpline(data['x'],data['y'][2],data['y'][3])
    def evaluate(x):return np.array([f(x),f(x,1),s(x),s(x,1)])
    record=dict(lambda_sigma=float(data['lambda_sigma']),beta=float(data['beta']),X=float(data['x'][-1]))
    return record,SimpleNamespace(sol=evaluate)

def main():
    out=HERE/'heldout';out.mkdir(exist_ok=False)
    prediction_bytes=(HERE/'predictions.json').read_bytes()
    predictions=json.loads(prediction_bytes)
    co=predictions['coefficients']
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),
        frozen_prediction_sha256=hashlib.sha256(prediction_bytes).hexdigest(),
        evaluations=[],tests=[])
    def save():(out/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    targets=predictions['declared_heldout_couplings']
    for kind in ['turn','binding']:
        seed=load_seed(kind)
        beta0=predictions['turn' if kind=='turn' else 'binding']['beta']
        for lam in [15.91,15.92,15.93,15.94,15.95,15.955,15.96,15.965,15.97,15.975]:
            r,sol=solve(lam,beta0,seed,X=40,tol=1e-9)
            report['evaluations'].append(r);save()
            if not r['accepted']:raise RuntimeError(str(r))
            seed=(r,sol)
            if lam not in targets:continue
            cache=[seed]
            def evaluate(beta):
                near=min(cache,key=lambda z:abs(z[0]['beta']-beta))
                rr,ss=solve(lam,float(beta),near,X=40,tol=1e-9)
                report['evaluations'].append(rr);save()
                if not rr['accepted']:raise RuntimeError(str(rr))
                cache.append((rr,ss));return rr,ss
            if kind=='turn':
                opt=minimize_scalar(lambda b:evaluate(b)[0]['Q'],bounds=(.123,.126),
                    method='bounded',options=dict(xatol=2e-7,maxiter=30))
                if not opt.success:raise RuntimeError(opt.message)
                r,sol=evaluate(float(opt.x))
            else:
                beta=brentq(lambda b:evaluate(b)[0]['binding_difference'],.264,.2665,xtol=1e-9)
                r,sol=evaluate(beta)
            delta=r['delta'];actual=r['scaled_charge']
            C0=co[kind+'_C0'];C1=co[kind+'_C1']
            pred=C0+delta*C1
            test=dict(kind=kind,lambda_sigma=lam,delta=delta,beta=r['beta'],Q=r['Q'],
                scaled_actual=actual,leading_prediction=C0,corrected_prediction=pred,
                leading_relative_error=abs(actual-C0)/actual,
                corrected_relative_error=abs(actual-pred)/actual,
                residual_over_delta_squared=(actual-pred)/(delta*delta),
                full_model_record=r)
            report['tests'].append(test);save()
            save_profile(out/(kind+'_lambda'+format(lam,'.6g')+'.npz'),r,sol)
            print(json.dumps(test),flush=True)
    report['source_sha256']={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ['heldout_test.py','full_solver.py']}
    save()
if __name__=='__main__':main()

