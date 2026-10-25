"""Build the held-out comparison table, checks and review figure."""
from pathlib import Path
import csv,hashlib,json,os
import numpy as np
HERE=Path(__file__).resolve().parent
os.environ['MPLCONFIGDIR']=str(HERE.parent.parent/'work'/'mpl-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
pred=json.loads((HERE/'predictions.json').read_text())
test=json.loads((HERE/'heldout/results.json').read_text())
validation=json.loads((HERE/'heldout_validation/results.json').read_text())
assert hashlib.sha256((HERE/'predictions.json').read_bytes()).hexdigest()==test['frozen_prediction_sha256']
for name,expected in pred['source_sha256'].items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==expected,name
rows=[]
for r in test['tests']:
    rows.append({k:r[k] for k in ['kind','lambda_sigma','delta','Q','scaled_actual','leading_prediction','corrected_prediction',
        'leading_relative_error','corrected_relative_error','residual_over_delta_squared']}|dict(
        error_improvement_factor=r['leading_relative_error']/r['corrected_relative_error']))
with (HERE/'heldout_comparison.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
summary=dict(status='Exploratory / Not Promoted',coefficients=pred['coefficients'],heldout_rows=rows,
    frozen_prediction_sha256=test['frozen_prediction_sha256'],
    max_validation_relative_charge_change=max(abs(r['relative_scaled_charge_change']) for r in validation['runs']),
    sensitivity_check=json.loads((HERE/'sensitivity_results.json').read_text())['turn_C1_difference_from_frozen_prediction'],
    publication_status='Working draft only; no submission or promotion')
(HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
fig,axes=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
for color,kind in [('C0','turn'),('C1','binding')]:
    rr=sorted([r for r in rows if r['kind']==kind],key=lambda r:r['delta'])
    x=[r['delta'] for r in rr]
    axes[0].loglog(x,[r['leading_relative_error'] for r in rr],'o--',color=color,label=kind+': leading')
    axes[0].loglog(x,[r['corrected_relative_error'] for r in rr],'s-',color=color,label=kind+': first correction')
    axes[1].plot(x,[r['residual_over_delta_squared'] for r in rr],'o-',color=color,label=kind)
axes[0].set(xlabel='delta',ylabel='Relative charge prediction error',title='Two declared held-out couplings')
axes[0].legend(fontsize=8)
axes[1].set(xlabel='delta',ylabel='(scaled charge - C0 - delta C1) / delta squared',title='Two-point remainder check')
axes[1].legend(fontsize=8)
fig.suptitle('Derived correction coefficients; no fit to full-model charge data',fontsize=12)
fig.savefig(HERE/'correction_evidence.png',dpi=160);plt.close(fig)
manifest={}
for f in sorted(HERE.rglob('*')):
    if f.is_file() and f.suffix in ['.py','.json','.csv','.npz','.png','.md','.txt'] and f.name!='manifest.json':
        manifest[str(f.relative_to(HERE)).replace('\\','/')]=hashlib.sha256(f.read_bytes()).hexdigest()
(HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(summary),flush=True)

