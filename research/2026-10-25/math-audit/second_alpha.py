"""Check the derived correction at a second alpha, without refitting it."""
import json
import numpy as np
from continue_math import ROOT,MC,central_solve,reduced_alpha,full_alpha
previous=None
for c in [.5,.7,.9,1.1,1.3]:
    previous,row=central_solve(c,320,1800,1e-9,previous)
    assert row['admissible']
a=float(previous.p[0])
red,base=reduced_alpha(a,previous)
h=5e-6
_,plus=reduced_alpha(a+h,red);_,minus=reduced_alpha(a-h,red)
q1=2*MC*(plus['S1']-minus['S1'])/(2*h)+(3/64-a)/MC*base['norm']/2
pred=dict(alpha=a,S0=base['S0'],S1=base['S1'],q0=base['q_scaled'],q1=q1)
(ROOT/'second_alpha_predictions_before_full.json').write_text(json.dumps(pred,indent=2))
rows=[]
for d in [.0625,.03125,.015625]:
    _,row=full_alpha(d,a,red)
    row.update(q_leading_error=row['q_scaled']-pred['q0'],
               q_first_order_error=row['q_scaled']-(pred['q0']+d*q1),
               S_leading_error=row['action_scaled']-pred['S0'],
               S_first_order_error=row['action_scaled']-(pred['S0']+d*pred['S1']))
    rows.append(row);print(json.dumps(row),flush=True)
assert all(abs(r['q_first_order_error'])<abs(r['q_leading_error']) for r in rows)
(ROOT/'second_alpha_results.json').write_text(json.dumps(dict(predictions=pred,full=rows),indent=2))
