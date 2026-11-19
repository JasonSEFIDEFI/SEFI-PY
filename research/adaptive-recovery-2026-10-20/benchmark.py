"""Frozen synthetic benchmark. Author project: Jason Duran Dutton."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import argparse,csv,hashlib,json,platform,sys,time,traceback
import numpy as np
from adaptive import predict,R
ROOT=Path(__file__).resolve().parent
P=json.loads((ROOT/'protocol.json').read_text())

def generate(seed,scenario):
    sid=P['scenarios'].index(scenario);streams=np.random.SeedSequence([seed,sid]).spawn(4)
    cr,br,mr,rr=[np.random.default_rng(x) for x in streams]
    n=P['stream_length'];u=np.linspace(0,1,n);mu=np.zeros((n,2));angle=np.zeros(n);sigma=np.full(n,.15)
    if scenario=='mean_step':
        mu[n//3:2*n//3]=(.45,-.30);mu[2*n//3:]=(-.25,.35)
    elif scenario=='rotation':angle=2*np.pi*u
    elif scenario=='joint_drift':
        mu=np.column_stack((.35*np.sin(4*np.pi*u),.25*np.cos(3*np.pi*u)-.25));angle=3*np.pi*u;sigma=.10+.18*(.5+.5*np.sin(2*np.pi*u))
    z=br.normal(size=(n,2))*np.column_stack((sigma,sigma*.5));c=np.cos(angle);s=np.sin(angle)
    beta=mu+np.column_stack((c*z[:,0]-s*z[:,1],s*z[:,0]+c*z[:,1]))
    if scenario=='bursts':beta+=(rr.random(n)<.05)[:,None]*rr.normal(0,.7,(n,2))
    y=beta+mr.normal(0,np.sqrt(R),(n,2))
    cal=cr.normal(size=(P['calibration_samples'],2))*[.15,.075]+cr.normal(0,np.sqrt(R),(P['calibration_samples'],2))
    return cal,y,beta

def fidelity(beta,estimate):
    r2=np.sum((beta-estimate)**2,axis=1);return np.exp(-r2)*(1-r2/2)**2
def savecsv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def interval(values):
    v=np.asarray(values);rng=np.random.default_rng(832719)
    boot=v[rng.integers(0,len(v),(10000,len(v)))].mean(axis=1)
    return dict(mean=float(v.mean()),lo=float(np.quantile(boot,.025)),hi=float(np.quantile(boot,.975)))

def run(stage):
    out=ROOT/'results'/stage;out.mkdir(parents=True,exist_ok=True)
    record={'status':'started','stage':stage,'protocol_sha256':digest(ROOT/'protocol.json'),'sources':{n:digest(ROOT/n) for n in ['adaptive.py','benchmark.py']},'python':sys.version,'numpy':np.__version__,'platform':platform.platform(),'processor':platform.processor(),'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
    manifest=out/'run_record.json';manifest.write_text(json.dumps(record,indent=2));started=time.perf_counter()
    try:
        if stage=='development':
            seeds=P['development_seeds'];methods=['ew:'+str(a) for a in P['baseline_grid']['ew_rates']]+['rolling:'+str(a) for a in P['baseline_grid']['rolling_windows']]
        else:
            lock=json.loads((ROOT/'selection.json').read_text());assert lock['protocol_sha256']==record['protocol_sha256']
            seeds=P['evaluation_seeds'];methods=['no_recovery','ml','fixed',lock['ew'],lock['rolling'],'mixture']
        rows=[];raw={};examples={};progress=[]
        for scenario in P['scenarios']:
            for seed in seeds:
                cal,y,beta=generate(seed,scenario)
                if stage=='evaluation':raw[f'{scenario}_{seed}_calibration']=cal;raw[f'{scenario}_{seed}_observations']=y;raw[f'{scenario}_{seed}_truth']=beta
                for method in methods:
                    ts=time.perf_counter();est=predict(method,cal,y);elapsed=time.perf_counter()-ts
                    assert est.shape==beta.shape and np.isfinite(est).all()
                    f=fidelity(beta,est);r2=np.sum((beta-est)**2,axis=1)
                    rows.append(dict(scenario=scenario,seed=seed,method=method,fidelity=float(f.mean()),residual_mse=float(r2.mean()),failure_fraction=float(np.mean(r2>.25**2)),seconds=elapsed))
                    if stage=='evaluation':raw[f'{scenario}_{seed}_{method}_estimate']=est
                    if stage=='evaluation' and seed==seeds[0]:examples[f'{scenario}_{method}']=f
                progress.append({'scenario':scenario,'seed':seed,'status':'completed'})
                record['completed_cases']=progress;manifest.write_text(json.dumps(record,indent=2))
            print(stage+' completed '+scenario,flush=True)
        savecsv(out/'seed_metrics.csv',rows)
        if stage=='development':
            scores={m:float(np.mean([r['fidelity'] for r in rows if r['method']==m])) for m in methods}
            ew=max((m for m in methods if m.startswith('ew:')),key=scores.get);rolling=max((m for m in methods if m.startswith('rolling:')),key=scores.get)
            selection={'protocol_sha256':record['protocol_sha256'],'ew':ew,'rolling':rolling,'primary_comparator':max([ew,rolling],key=scores.get),'development_scores':scores,'evaluation_not_used':True,'source_hashes':record['sources']}
            (ROOT/'selection.json').write_text(json.dumps(selection,indent=2));print(json.dumps(selection),flush=True)
        else:
            np.savez_compressed(out/'raw_streams_and_predictions.npz',**raw)
            np.savez_compressed(out/'example_fidelities.npz',**examples)
            summary=[];comparisons=[]
            for s in P['scenarios']:
                for m in methods:
                    vals=[r['fidelity'] for r in rows if r['scenario']==s and r['method']==m]
                    summary.append(dict(scenario=s,method=m,**interval(vals)))
                    if m!='mixture':
                        diffs=[next(r['fidelity'] for r in rows if r['scenario']==s and r['method']=='mixture' and r['seed']==k)-next(r['fidelity'] for r in rows if r['scenario']==s and r['method']==m and r['seed']==k) for k in seeds]
                        comparisons.append(dict(scenario=s,baseline=m,**interval(diffs)))
            savecsv(out/'fidelity_summary.csv',summary);savecsv(out/'paired_differences.csv',comparisons)
            primary=next(r for r in comparisons if r['scenario']==P['primary_scenario'] and r['baseline']==lock['primary_comparator'])
            decision={'primary':primary,'improvement_supported':primary['lo']>0,'candidate':'conventional adaptive mixture, no established SEFI-specific novelty','evaluation_cases':len(progress),'evaluation_samples':len(progress)*P['stream_length']}
            (out/'decision.json').write_text(json.dumps(decision,indent=2));print(json.dumps(decision),flush=True)
        record.update(status='completed',elapsed_seconds=time.perf_counter()-started,methods=methods)
    except BaseException as e:
        record.update(status='interrupted' if isinstance(e,KeyboardInterrupt) else 'failed',error=type(e).__name__)
        raise
    finally:manifest.write_text(json.dumps(record,indent=2))

def timing():
    lock=json.loads((ROOT/'selection.json').read_text());cal,y,_=generate(2001,'joint_drift')
    methods=['fixed',lock['ew'],lock['rolling'],'mixture'];rows=[]
    for m in methods:predict(m,cal,y[:200])
    for rep in range(7):
        for m in (methods if rep%2==0 else list(reversed(methods))):
            t=time.perf_counter();predict(m,cal,y);elapsed=time.perf_counter()-t
            rows.append(dict(repetition=rep,method=m,seconds=elapsed,microseconds_per_sample=elapsed/len(y)*1e6))
    out=ROOT/'results/evaluation';savecsv(out/'timing.csv',rows)
    print(json.dumps({m:float(np.median([r['microseconds_per_sample'] for r in rows if r['method']==m])) for m in methods}),flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['development','evaluation','timing']);args=a.parse_args()
    timing() if args.stage=='timing' else run(args.stage)
