from pathlib import Path
import json,hashlib,subprocess
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parent
repo=ROOT.parents[2]
record=dict(status='stopped_at_user_request',recorded_utc=datetime.now(timezone.utc).isoformat(),
    authorized_limit_minutes=30,stop_reason='User requested stopping early and recording the vortex work; compute was not exhausted.',
    completed_cases=['pilot_control_n8','pilot_m2_n8','nonlinear_m2_n16'],
    interrupted_case='nonlinear_m2_n24',interrupted_case_state_saved=False,
    not_run=['nonlinear_m2_n16_half_dt','nonlinear_strong_m2_n16','nonlinear_m2_wide_n16','long_linear'],
    limits=['No angular-resolution convergence comparison completed','No temporal-refinement comparison completed for nonlinear angular extension','No independent Authorship coordinate evolved','No proof of isolated matter or nonlinear stability'],
    validation={'prototype_consistency_tests':4,'new_evidence_route_test':1,'dispersion_route_math_tests':5},
    protected_equations={})
for name in ['physics.py','ring_operator.py','rest_search.py']:
    file='SEFI-PY Evolution/sim/vortex_console/'+name
    original=subprocess.check_output(['git','-C',str(repo),'show','HEAD:'+file]).replace(b'\r\n',b'\n')
    current=(repo/file).read_bytes().replace(b'\r\n',b'\n')
    assert current==original
    record['protected_equations'][name]=hashlib.sha256(current).hexdigest()
(ROOT/'STOP_RECORD.json').write_text(json.dumps(record,indent=2)+'\n')
(ROOT/'.gitignore').write_text('__pycache__/\n*.pyc\n')
(ROOT/'README.md').write_text('''# Experimental nonlinear angular extension - paused checkpoint

Jason Duran Dutton | 19 October 2026

The user stopped this batch early. It did not exhaust the authorized 30-minute allowance. See STOP_RECORD.json for complete, interrupted and unrun cases. No simulation is scheduled to resume automatically.

The separate cylindrical prototype extends the preserved lambda=25, N=0 equations into an angularly periodic, nonlinear r/theta/Z evolution. It uses finite differences in r,Z and a Fourier angular Laplacian, with axis regularity enforced by a single angularly averaged axis value. It is not an Authorship-dependent model. The main console's physics.py, ring_operator.py and rest_search.py are unchanged.

Four prototype tests pass: reduction to the axisymmetric equations; a harmonic r^2 exp(2i theta) angular-Laplacian check; directional agreement with the prior angular tangent equations; and differential discrete energy/carrier-charge checks. These are consistency checks, not empirical validation.

Completed pilot control: 8 angles, T=8, no imposed disturbance; all sampled meridional angles retain a detected core. Pilot m=2 disturbance: 8 angles, T=8. Main m=2 disturbance: 16 angles, T=12; all sampled angles retain a detected core, maximum angular asymmetry gain 1.2084, relative energy drift 3.23e-11 and carrier-charge drift 2.77e-11.

The 24-angle run was interrupted without saving its in-progress state. Angular and temporal convergence have therefore NOT been established. Pilot and main runs have different final times and cannot substitute for a matched refinement comparison. The selected m=2 initial condition retains a discrete symmetry and does not test arbitrary 3D disturbances. Nonlinear collocation products can alias; more angles, spectral diagnostics and matched runs are needed. Fixed boundaries, finite time, incomplete spatial convergence and absence of Authorship coupling remain substantive limitations.

Report.pdf records both this partial batch and the earlier completed survival study. The console receives a read-only evidence snapshot; experimental solvers do not become selectable production models.

Reproduction: install NumPy; run `python -m unittest test_extended -v`. Invoke `cylindrical3d.run_case` with explicit parameters for an individual case. `batch.py` preserves the dated plan and original cutoff; do not assume it is a currently scheduled computation. `snapshot.py` regenerates the read-only evidence snapshot from saved completed cases. `record_stop.py` verifies protected equations and records this stop checkpoint. Pillow and ReportLab are needed only for figures/reporting.
''',encoding='utf-8')
print('Stop checkpoint recorded')
