"""Independent domain/tolerance refinement of the two reduced landmarks."""
import json
from scipy.optimize import minimize_scalar,brentq
from continue_math import ROOT,central_solve,integrate_profile

previous=None
for c in [.5,.7,.9,1.1,1.3,1.45]:
    previous,row=central_solve(c,320,1800,1e-9,previous)
    assert row['admissible']
seed=previous
cache={}
def at(c):
    key=round(float(c),13)
    if key not in cache:
        sol,row=central_solve(float(c),400,2600,1e-10,seed)
        assert row['admissible'],row
        cache[key]=integrate_profile(sol,row['alpha'],400)
    return cache[key]
opt=minimize_scalar(lambda c:at(c)['q_scaled'],bounds=(1.06,1.15),method='bounded',options={'xatol':1e-8})
root=brentq(lambda c:at(c)['binding_scaled'],1.4,1.44,xtol=1e-11)
result=dict(L=400,tol=1e-10,fold=at(opt.x),energetic_crossing=at(root))
(ROOT/'landmark_refinement.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
