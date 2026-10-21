"""Independent validation of saved outcomes; no rerunning/tuning estimators."""
from pathlib import Path
import csv,hashlib,json
import numpy as np
from physical_reference import metrics
P=Path(__file__).resolve().parent
protocol=json.loads((P/'protocol.json').read_text());selection=json.loads((P/'selection.json').read_text())
record=json.loads((P/'results/evaluation/run_record.json').read_text())
assert record['status']=='completed'
assert record['protocol_sha256']==hashlib.sha256((P/'protocol.json').read_bytes()).hexdigest()==selection['protocol_sha256']
for name,sha in record['sources'].items():assert hashlib.sha256((P/name).read_bytes()).hexdigest()==sha==selection['source_hashes'][name]
rows=list(csv.DictReader((P/'results/evaluation/seed_metrics.csv').open()))
data=np.load(P/'results/evaluation/raw_streams_and_predictions.npz');maximum=0.
assert len(rows)==len(protocol['evaluation_seeds'])*len(protocol['scenarios'])*6
for row in rows:
    key=row['scenario']+'_'+row['seed'];beta=data[key+'_truth'];estimate=data[key+'_'+row['method']+'_estimate']
    independent=metrics(beta-estimate)['entanglement_fidelity'].mean()
    maximum=max(maximum,abs(float(independent)-float(row['fidelity'])))
assert maximum<1e-13
result={'status':'passed','rows_verified':len(rows),'complete_cases':len(record['completed_cases']),'scored_samples_including_six_methods':len(rows)*protocol['stream_length'],'max_fidelity_discrepancy':maximum,'source_and_protocol_hashes_match_development_lock':True,'incomplete_cases':sum(x['status']!='completed' for x in record['completed_cases'])}
(P/'results/evaluation/saved_results_validation.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
