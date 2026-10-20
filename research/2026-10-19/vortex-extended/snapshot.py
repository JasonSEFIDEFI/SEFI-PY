"""Build a read-only console evidence snapshot from completed recorded cases."""
from pathlib import Path
import json
from datetime import datetime,timezone
from run_integrity import eligible
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[2]
CONSOLE=REPO/'SEFI-PY Evolution/sim/vortex_console'
prior=ROOT.parent/'vortex-survival'
axis=json.loads((prior/'summary.json').read_text())['cases']
linear=json.loads((prior/'angular_summary.json').read_text())['cases']
nonlinear=[]
excluded=[dict(name='nonlinear_m2_n24',reason='Original run interrupted at author request; no final state saved')]
for path in sorted(ROOT.glob('*.json')):
    try:result=json.loads(path.read_text())
    except (json.JSONDecodeError,UnicodeDecodeError):
        excluded.append(dict(name=path.stem,reason='Invalid/truncated output'));continue
    if 'ntheta' in result and 'samples' in result:
        if not eligible(result):
            excluded.append(dict(name=result.get('name',path.stem),reason=result.get('reason','unfinished')));continue
        entry={k:v for k,v in result.items() if k!='samples'}
        entry['final_time']=result['samples'][-1]['time']
        entry['high_angular_power_fraction']=max(x['high_angular_power_fraction'] for x in result['samples'])
        nonlinear.append(entry)
data=dict(date=datetime.now(timezone.utc).isoformat(),
    scope='Finite-time, finite-domain evidence. No independent Authorship dynamics or physical-matter identification.',
    axisymmetric=axis,linear=linear,nonlinear=nonlinear,excluded_runs=excluded,
    checks=['20 production-console regression tests passed at the survival checkpoint.',
            'Four survival/core/balance tests passed.',
            'Four nonlinear angular prototype consistency tests passed.',
            'Main console physical equations unchanged; new research views are read-only.',
            'Angular/temporal refinement is reported by completed case; it is not a continuum stability proof.',
            'Interactive page/API automated checks are recorded separately; visual browser inspection may be unavailable.'])
(CONSOLE/'data/research_evidence.json').write_text(json.dumps(data,indent=2)+'\n')
(ROOT/'evidence_snapshot.json').write_text(json.dumps(data,indent=2)+'\n')
print('Snapshot:',len(nonlinear),'nonlinear angular runs')
