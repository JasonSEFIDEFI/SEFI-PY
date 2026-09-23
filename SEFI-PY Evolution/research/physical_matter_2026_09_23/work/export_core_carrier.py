import json
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
names=['core_condensate_L20_n400','core_condensate_L20_n800','core_condensate_L30_n600','core_condensate_L20_n400_omega0.3','core_condensate_L20_n800_omega0.3']
records=[json.loads((root/'work'/f'{name}.json').read_text()) for name in names]
report={'status':'Provisional two-field classical model; no physical-particle identification or full stability proof.', 'parameters':{'alpha':1,'beta':1,'Omega_phi':float(np.sqrt(2)),'kappa':4,'mu_sigma_squared':3,'lambda_sigma':25},'profiles':records}
(root/'outputs'/'core_carrier_nonlinear_results.json').write_text(json.dumps(report,indent=2))
p=np.load(root/'work'/'core_condensate_L20_n800_omega0.3.npz')
np.savetxt(root/'outputs'/'core_carrier_charged_profile.csv',np.column_stack([p['r'],p['f'],p['s']]),delimiter=',',header='radius,original_field_amplitude,carrier_amplitude',comments='')
print(json.dumps({'charged_refined':records[-1],'relative_norm_refinement':abs(records[-1]['carrier_transverse_norm']/records[-2]['carrier_transverse_norm']-1)},indent=2))
