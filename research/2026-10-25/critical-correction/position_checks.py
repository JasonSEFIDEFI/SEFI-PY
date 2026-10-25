"""Secondary reduced-equation checks of beta-root shifts, after charge tests.
These position coefficients were not part of the frozen held-out protocol.
"""
from pathlib import Path
import json
import numpy as np
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    out=HERE/'position_checks.json'
    if out.exists():raise FileExistsError(out)
    p=json.loads((HERE/'predictions.json').read_text())
    s=json.loads((HERE/'sensitivity_results.json').read_text())
    h=json.loads((HERE/'heldout/results.json').read_text())
    mc=np.sqrt(1.75);b=p['turn']['beta'];center=p['turn']
    checks=[]
    for r in p['derivative_checks']:
        a,lo,hi,d=r['evaluations'];step=r['step']
        second=lambda key:(-a[key]+16*lo[key]-30*center[key]+16*hi[key]-d[key])/(12*step*step)
        I0pp=second('I0');A1pp=second('A1')
        Q0pp=mc*I0pp
        I0p=s['runs'][0]['I0_beta_derivative']
        Q1p=2*mc*A1pp+((2-b)*I0p-center['I0'])/(2*mc)
        checks.append(dict(step=step,Q0_second_derivative=Q0pp,Q1_beta_derivative=Q1p,
            beta_turn_first_correction=-Q1p/Q0pp))
    beta_turn1=checks[-1]['beta_turn_first_correction']
    beta_binding1=(p['coefficients']['binding_C1']-s['runs'][1]['Q1'])/s['runs'][1]['Q0_beta_derivative']
    tests=[]
    for test in h['tests']:
        kind=test['kind'];base=p['turn' if kind=='turn' else 'binding']['beta']
        correction=beta_turn1 if kind=='turn' else beta_binding1
        predicted=base+test['delta']*correction
        tests.append(dict(kind=kind,lambda_sigma=test['lambda_sigma'],actual_beta=test['beta'],
            first_order_beta=predicted,
            remainder_over_delta_squared=(test['beta']-predicted)/test['delta']**2))
    result=dict(scope=__doc__,turn_step_checks=checks,beta_turn_C1=beta_turn1,
        beta_binding_C1=beta_binding1,tests=tests)
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

