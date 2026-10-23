"""Fresh finite-grid fixed-charge relaxation and curvature runs of existing model."""
from pathlib import Path
import json
import sys
import time
ROOT=Path(__file__).resolve().parents[3]
SOURCE=ROOT/'research/2026-10-19/vortex-extended'
sys.path.insert(0,str(SOURCE))
import vacuum_spherical as model
import spherical_stability as stability

CONFIG=dict(grids=[dict(h=.8,L=60),dict(h=.4,L=60),dict(h=.2,L=60),dict(h=.2,L=80)],
            wall_seconds_per_relaxation=20,maximum_iterations=20000,
            note='Uncalibrated dimensionless two-field model, Q=100000, lambda=12. No empirical particle identification.')

def main():
    out=Path(__file__).parent/'results';out.mkdir(exist_ok=True)
    model.ROOT=out;stability.ROOT=out
    (out/'config.json').write_text(json.dumps(CONFIG,indent=2))
    errors=[]
    for grid in CONFIG['grids']:
        tag=f"h{grid['h']}_L{grid['L']}"
        try:
            model.run(**grid,wall=CONFIG['wall_seconds_per_relaxation'],maximum=CONFIG['maximum_iterations'],tag=tag)
            meta=json.loads((out/(tag+'.json')).read_text())
            if meta['converged']: stability.run(tag)
            else: errors.append(dict(tag=tag,error='Not converged; curvature calculation skipped'))
        except Exception as exc: errors.append(dict(tag=tag,error=repr(exc)))
        (out/'failures.json').write_text(json.dumps(errors,indent=2))
    print(json.dumps(dict(failures=errors)))

if __name__=='__main__':main()
