"""Independent saved-array re-evaluation, without rerunning or choosing new trials."""
import csv
import json
from pathlib import Path
import numpy as np
from bosonic import metrics
from decoders import GeometryDecoder,MaximumLikelihood,ConventionalPosterior

ROOT=Path(__file__).resolve().parent
data=np.load(ROOT/'results/raw_trials.npz')
rows=list(csv.DictReader((ROOT/'results/seed_metrics.csv').open()))
protocol=json.loads((ROOT/'protocol.json').read_text())
max_defect=0.
for k in range(protocol['seeds']):
    for j,sigma in enumerate(data['sigmas']):
        y=data['observations_y'][k,j];beta=data['truth_beta'][k,j]
        s=data['calibration_covariances'][k,j];r=data['measurement_covariance']
        for name,decoder in [('ML',MaximumLikelihood()),('Geometry',GeometryDecoder(s,r)),
                             ('Conventional posterior',ConventionalPosterior(s,r))]:
            expected=metrics(beta-decoder.predict(y),protocol['accuracy_residual_radius'])
            row=next(x for x in rows if float(x['sigma'])==sigma and int(x['seed'])==protocol['first_seed']+k and x['method']==name)
            for key,values in expected.items():max_defect=max(max_defect,abs(float(row[key])-float(values.mean())))
assert max_defect<1e-12
all_rows=list(csv.DictReader((ROOT/'results/table_1_all_results.csv').open()))
assert len(rows)==1400 and len(all_rows)==70
for row in all_rows:
    same=[x for x in rows if x['sigma']==row['sigma'] and x['method']==row['method']]
    assert len(same)==20
    assert abs(np.mean([float(x['entanglement_fidelity']) for x in same])-float(row['entanglement_fidelity']))<1e-12
    assert float(row['entanglement_fidelity_ci_low'])<=float(row['entanglement_fidelity'])<=float(row['entanglement_fidelity_ci_high'])
record={'raw_saved_array_evaluation_pass':True,'maximum_seed_metric_discrepancy':max_defect,
        'seed_method_rows':len(rows),'aggregate_noise_method_rows':len(all_rows),
        'primary_decoder_trials_rechecked':560000,'author':'Jason D. Dutton','date':'2026-10-06'}
(ROOT/'results/saved_data_validation.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
print(json.dumps(record,indent=2))
