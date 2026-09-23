"""Algebraic audit of the time-harmonic rest obstruction and an escape window."""
import json
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
rng=np.random.default_rng(1729);x=rng.uniform(0,8,10000);y=rng.uniform(0,4,10000);nu=rng.uniform(0,1,10000)
direct=(x-1)**2/4+25*y*y/4+(4*x-3-nu*nu)*y/2
squares=((x-1)+4*y)**2/4+9*y*y/4+(1-nu*nu)*y/2
lam=12.;trial_nu=.85;yy=(3+trial_nu**2)/lam
result={'identity_max_abs_error':float(abs(direct-squares).max()),'baseline_lambda':25,'baseline_result':'For |nu|<=1, U_nu is nonnegative. Derrick scaling K+3U=0 excludes nontrivial finite-relative-energy, localized time-harmonic rest profiles on R^3 with uniform background and canonical gradients. Does not exclude moving loops, extended strings or more general time dependence.', 'proposed_parameter_change':{'lambda':lam,'nu':trial_nu,'bulk_mass_squared':1,'zero_frequency_global_vacuum_condition':'lambda >= 9 (here kappa=4, mu^2=3, beta=1, R_infinity=1)','open_escape_window':'9 < lambda < 16 and sqrt(lambda)-3 < nu^2 < 1','negative_effective_potential_trial':{'phi_squared':0,'sigma_squared':yy,'U_nu':float(.25-(3+trial_nu**2)**2/(4*lam))},'interpretation':'Removes this positivity obstruction only. Does not prove existence, stability, or physical identity.'}}
assert result['identity_max_abs_error']<1e-12
(root/'outputs'/'rest_branch_audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
