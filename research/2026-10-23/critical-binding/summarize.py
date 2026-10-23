import json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from critical import ROOT

data=json.loads((ROOT/'validation.json').read_text())
rows=data['full'];ref=data['reference']
assert len(rows)==8 and all(r['success'] for r in rows)
errors=np.abs([r['relative_q_error'] for r in rows])
assert np.all(np.diff(errors)<0)
fig,axes=plt.subplots(1,2,figsize=(10,4))
d=np.array([r['delta'] for r in rows]);q=np.array([r['q_scaled'] for r in rows])
axes[0].semilogx(d,q,'o-',label='Original two-field equations')
axes[0].axhline(ref['q_scaled'],color='k',ls='--',label='Reduced prediction')
axes[0].set(xlabel='16 - lambda',ylabel='Q (16 - lambda)^(3/2)',title='Fixed central amplitude; not a minimum-charge fit')
axes[0].legend(fontsize=8)
axes[1].semilogx(d,[r['binding_scaled'] for r in rows],'o-')
axes[1].axhline(ref['binding_scaled'],color='k',ls='--')
axes[1].axhline(0,color='gray',lw=.7)
axes[1].set(xlabel='16 - lambda',ylabel='(E - m Q) sqrt(16 - lambda)',title='Energetic dispersion threshold at zero')
fig.tight_layout();fig.savefig(ROOT/'critical_convergence.png',dpi=180);plt.close(fig)
last=rows[-1]
report=f'''# Observed outcomes

Exact symbolic identities passed. The initial scan accepted nine of ten solves; the alpha<0 case is excluded. Six domain/tolerance checks and eight original two-field solves converged. A failed direct-initialization attempt is retained separately.

At central angle 1.45 the reduced prediction is scaled charge {ref['q_scaled']:.8f}, scaled binding {ref['binding_scaled']:.8f}. The full-model charge discrepancy decreases monotonically from {100*errors[0]:.3f}% at delta=1 to {100*errors[-1]:.3f}% at delta={last['delta']}. This supports the formal asymptotic reduction on this branch. It does not measure the minimum-charge scaling exponent independently.

At this central amplitude the full-model scaled binding changes sign between delta=0.25 and 0.125. Thus binding cannot be inferred quantitatively from the leading critical limit at delta=1. No Hessian or dynamical stability result is supplied in this experiment.

For central angle 1.45, the L=240,tolerance=1e-8 and L=320,tolerance=1e-9 reduced scaled charges differ by less than 1e-8 relatively. Full-model domain refinement remains to be done. Residuals, timings, configurations and profiles are saved alongside this report. This is deterministic numerical convergence, not a statistical uncertainty estimate.

No experimental discovery, new element, unique SEFI advantage, or priority claim follows. See MATH.md for the derivation, assumptions and next tests.
'''
(ROOT/'OUTCOMES.md').write_text(report)
files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.iterdir() if p.is_file() and p.name!='SHA256SUMS.json'}
(ROOT/'SHA256SUMS.json').write_text(json.dumps(files,indent=2))
print(report)
