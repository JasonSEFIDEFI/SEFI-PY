from pathlib import Path
import json,csv,hashlib,math
from phase_model import simulate,equilibrium,critical_detuning
P=Path(__file__).resolve().parent
def csvwrite(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def run():
    protocol=json.loads((P/'protocol.json').read_text());out=P/'results';out.mkdir(exist_ok=True)
    record={'status':'started','protocol_sha256':hashlib.sha256((P/'protocol.json').read_bytes()).hexdigest(),'source_sha256':hashlib.sha256((P/'phase_model.py').read_bytes()).hexdigest(),'completed_cases':[]}
    manifest=out/'run_record.json';manifest.write_text(json.dumps(record,indent=2));rows=[];traces=[]
    try:
        for name,alpha in [('constant_gain',0.),('gain_degradation',protocol['alpha_test'])]:
            for delta in protocol['detunings']:
                for dt in protocol['dt_grid']:
                    r,trace=simulate(delta,protocol['K0'],alpha,protocol['beta'],dt,protocol['duration_seconds'],protocol['late_window_seconds'])
                    eq=equilibrium(delta,protocol['K0'],alpha,protocol['beta'])
                    rows.append(dict(model=name,detuning=delta,dt=dt,alpha=alpha,**r,equilibrium_exists=eq is not None,predicted_phase=eq[0] if eq else '',predicted_degradation=eq[1] if eq else ''))
                    assert r['min_degradation']>-1e-12
                    if eq:
                        assert abs(r['final_phase']-eq[0])<1e-7 and abs(r['final_degradation']-eq[1])<1e-7
                    if dt==min(protocol['dt_grid']):
                        for t,p,d in trace:traces.append(dict(model=name,detuning=delta,t=t,phase=p,degradation=d))
                    record['completed_cases'].append([name,delta,dt]);manifest.write_text(json.dumps(record,indent=2))
        csvwrite(out/'all_trajectories_summary.csv',rows);csvwrite(out/'sampled_trajectories.csv',traces)
        thresholds=[dict(alpha=a,critical_detuning=critical_detuning(protocol['K0'],a,protocol['beta'])) for a in protocol['alpha_threshold_grid']]
        csvwrite(out/'locking_thresholds.csv',thresholds)
        convergence=[]
        for name in ['constant_gain','gain_degradation']:
            for delta in protocol['detunings']:
                group=[r for r in rows if r['model']==name and r['detuning']==delta];fine=min(group,key=lambda r:r['dt'])
                convergence.append(dict(model=name,detuning=delta,max_final_phase_difference=max(abs(r['final_phase']-fine['final_phase']) for r in group),all_lock_diagnostics_agree=len({r['numerically_locked'] for r in group})==1))
        csvwrite(out/'time_step_convergence.csv',convergence)
        record.update(status='completed',cases=len(rows),lock_classifications_agree_all_steps=all(x['all_lock_diagnostics_agree'] for x in convergence),max_final_phase_difference=max(x['max_final_phase_difference'] for x in convergence))
        print(json.dumps(record,indent=2));print(json.dumps(thresholds,indent=2))
    except BaseException as e:
        record.update(status='interrupted' if isinstance(e,KeyboardInterrupt) else 'failed',error=type(e).__name__);raise
    finally:manifest.write_text(json.dumps(record,indent=2))
if __name__=='__main__':run()
