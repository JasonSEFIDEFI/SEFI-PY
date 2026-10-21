"""Exploratory physical-fidelity intervention rule; conventional Bayesian decision theory."""
from pathlib import Path
import numpy as np
import json,csv,time,hashlib
P=Path(__file__).resolve().parent

def expected_fidelity(mean,cov):
    """Exact Gaussian average of exp(-r²)(1-r²/2)²; derived in DERIVATION.md."""
    mean=np.atleast_2d(mean);v,u=np.linalg.eigh(cov);assert v.min()>-1e-12
    z=mean@u;inv=1/(1+2*v)
    g=np.prod(inv)**.5*np.exp(-np.sum(z*z*inv,axis=1))
    l1=-np.sum(v*inv)-np.sum(z*z*inv**2,axis=1)
    l2=2*np.sum(v*v*inv**2)+4*np.sum(z*z*v*inv**3,axis=1)
    return g*(1+l1+.25*(l1*l1+l2))

def score(residual):
    r2=np.sum(residual**2,axis=1);return np.exp(-r2)*(1-r2/2)**2

def independent_quadrature(mean,cov):
    x,w=np.polynomial.hermite.hermgauss(40);i,j=np.meshgrid(np.arange(40),np.arange(40))
    z=np.column_stack((x[i.ravel()],x[j.ravel()]))*np.sqrt(2)
    return np.sum(score(np.asarray(mean)+z@np.linalg.cholesky(cov).T)*(w[i.ravel()]*w[j.ravel()])/np.pi)

def run():
    p=json.loads((P/'protocol.json').read_text());out=P/'results';out.mkdir(exist_ok=True)
    record={'status':'started','protocol_sha256':hashlib.sha256((P/'protocol.json').read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'completed_seeds':[]}
    (out/'run_record.json').write_text(json.dumps(record,indent=2))
    try:
        # Numerical integration checks a full correlated covariance, not only diagonal examples.
        check_cov=np.array([[.04,.01],[.01,.02]]);error=0.
        for mean in [[0,0],[.3,-.2],[1,.5]]:
            error=max(error,abs(float(expected_fidelity(mean,check_cov)[0])-independent_quadrature(mean,check_cov)))
        assert error<1e-11
        s=np.array(p['prior_covariance']);r=np.array(p['measurement_covariance']);c=np.linalg.inv(np.linalg.inv(s)+np.linalg.inv(r));w=s@np.linalg.inv(s+r)
        rows=[];raw={};max_risk_violation=0.
        for seed in range(p['seeds']):
            rng=np.random.default_rng(81000+seed);n=p['trials_per_seed'];beta=rng.multivariate_normal([0,0],s,n);y=beta+rng.multivariate_normal([0,0],r,n);mu=y@w.T
            raw[f'{seed}_beta']=beta;raw[f'{seed}_observations']=y
            hold_score=score(beta);hold_expected=expected_fidelity(mu,c)
            for q in p['actuator_variance_per_real_component']:
                disturbance=rng.normal(0,np.sqrt(q),(n,2));act_score=score(beta-mu+disturbance)
                act_expected=float(expected_fidelity([0,0],c+q*np.eye(2))[0])
                gate=act_expected>hold_expected;mse=np.sum(mu*mu,axis=1)>2*q
                max_risk_violation=max(max_risk_violation,float(np.max(np.maximum(act_expected,hold_expected)-np.where(gate,act_expected,hold_expected))))
                raw[f'{seed}_{q}_actuator_noise']=disturbance
                for name,chosen in [('hold',np.zeros(n,dtype=bool)),('always_correct',np.ones(n,dtype=bool)),('mse_gate',mse),('fidelity_gate',gate)]:
                    achieved=np.where(chosen,act_score,hold_score)
                    predicted=np.where(chosen,act_expected,hold_expected)
                    rows.append({'seed':seed,'actuator_variance':q,'policy':name,'fidelity':float(achieved.mean()),'predicted_fidelity':float(predicted.mean()),'action_fraction':float(chosen.mean())})
            record['completed_seeds'].append(seed);(out/'run_record.json').write_text(json.dumps(record,indent=2))
        with (out/'seed_results.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        np.savez_compressed(out/'simulation_inputs.npz',**raw)
        summary=[]
        for q in p['actuator_variance_per_real_component']:
            values={name:np.array([z['fidelity'] for z in rows if z['actuator_variance']==q and z['policy']==name]) for name in p['policies']}
            entry={'actuator_variance':q,'policies':{name:{'fidelity':float(v.mean()),'action_fraction':float(np.mean([z['action_fraction'] for z in rows if z['actuator_variance']==q and z['policy']==name]))} for name,v in values.items()},'paired_gain':{}}
            for baseline in ['always_correct','mse_gate','hold']:
                delta=values['fidelity_gate']-values[baseline];boot=np.random.default_rng(883).choice(delta,(10000,len(delta)),replace=True).mean(axis=1)
                entry['paired_gain'][baseline]={'mean':float(delta.mean()),'lo':float(np.quantile(boot,.025)),'hi':float(np.quantile(boot,.975))}
            summary.append(entry)
        (out/'summary.json').write_text(json.dumps(summary,indent=2))
        record.update(status='completed',quadrature_max_error=error,conditional_decision_violation=max_risk_violation,unique_displacement_trials=p['seeds']*n,all_four_actuator_settings_reported=True)
        print(json.dumps(summary,indent=2),flush=True)
    except BaseException as exc:
        record.update(status='interrupted' if isinstance(exc,KeyboardInterrupt) else 'failed',error=type(exc).__name__);raise
    finally:(out/'run_record.json').write_text(json.dumps(record,indent=2))

if __name__=='__main__':run()
