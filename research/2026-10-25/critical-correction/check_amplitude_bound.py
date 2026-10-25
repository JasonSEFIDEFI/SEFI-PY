"""Exact algebra checks for the conditional amplitude estimate."""
from pathlib import Path
import json
import sympy as sp
HERE=Path(__file__).resolve().parent
d,lam,s0,chi,lap0,lapchi,c=sp.symbols('d lam s0 chi lap0 lapchi c',real=True)
sapp=s0+d*chi
raw=-d*lap0-d*d*lapchi+lam*(sapp*sapp-s0*s0)*sapp
expanded=d*(-lap0+2*lam*s0*s0*chi)+d*d*(-lapchi+3*lam*s0*chi*chi)+d**3*lam*chi**3
identity=sp.expand(raw-expanded)
residual=sp.expand((raw-d*d*(-lapchi+3*lam*s0*chi*chi)-d**3*lam*chi**3).subs(lap0,2*lam*s0*s0*chi))
coercivity=sp.simplify((3*sp.Rational(3,4)**2-1)*3*c*c-sp.Rational(33,16)*c*c)
assert identity==0 and residual==0 and coercivity==0
result=dict(scope=__doc__,residual_identity=str(identity),leading_cancellation=str(residual),
    coercivity_coefficient=str(coercivity),
    limitation='Symbolic identities only; no certified verification of branch hypotheses.')
out=HERE/'amplitude_bound_checks.json'
if out.exists():raise FileExistsError(out)
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))

