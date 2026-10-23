"""Locked extension: stronger residual models and calibration under chart shift.

All tests are retrospective; no claim of independent preregistration is made.
The protocol is written before fitting. Test targets never select hyperparameters.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import platform
import sys
import time
import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.spatial.distance import cdist

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('pilot',ROOT/'analyze.py')
pilot=importlib.util.module_from_spec(spec);spec.loader.exec_module(pilot)
spec=importlib.util.spec_from_file_location('old_data',ROOT/'historical.py')
old_data=importlib.util.module_from_spec(spec);spec.loader.exec_module(old_data)
CONFIG=dict(version='publication-extension-1',random_seeds=[8400,8401,8402,8403,8404],
    historical_seeds=[9600,9601,9602,9603,9604],element_chain_folds=5,
    length_scales=[1.5,3.,6.],noise_fractions=[.03,.3],parity_coordinate_scale=2.,
    interval_alpha=.1,bootstrap_seed=2307,bootstrap_repeats=2000,
    primary_endpoint='Total binding-energy RMSE (MeV); actual held-out AME2016-to-2020 cohort',
    secondary_endpoints=['MAE','nominal-90% residual interval coverage','distance-scaled coverage',
        'Gaussian process predictive-normal and residual-calibrated coverage','width','setup and inference time'],
    model_families=['SEMF','SEMF_known_shells','Gaussian_residual_smoother',
        'GP_RBF','GP_Matern32','GP_RBF_parity'],
    interpretation='Retrospective extension following pilot; no SEFI-specific law or discovery claim.')

def load_arrays():
    _,current,audit=pilot.load()
    old=old_data.load_2016()
    older=[r for r in old.values() if r['Z']>0 and r['N']>0 and r['A']>=16
           and not r['estimated'] and r['binding_per_A_keV'] is not None]
    x=np.array([[r['Z'],r['N']] for r in current],float)
    y=np.array([r['binding_per_A_keV']*r['A']/1000 for r in current])
    xo=np.array([[r['Z'],r['N']] for r in older],float)
    yo=np.array([r['binding_per_A_keV']*r['A']/1000 for r in older])
    later=np.array([i for i,r in enumerate(current) if (r['Z'],r['N']) not in old
                    or old[(r['Z'],r['N'])]['estimated']],int)
    return x,y,xo,yo,later,audit

def partitions():
    x,y,xo,yo,later,audit=load_arrays();n=len(x)
    for seed in CONFIG['random_seeds']:
        order=np.random.default_rng(seed).permutation(n)
        a,b,c=np.split(order,[int(.6*n),int(.8*n)])
        va=np.random.default_rng(seed+1).permutation(b);tu,ca=np.array_split(va,2)
        yield f'random_{seed}',x[a],y[a],x[tu],y[tu],x[ca],y[ca],x[c],y[c]
    z=x[:,0]
    for fold in range(CONFIG['element_chain_folds']):
        tr=np.where(z%5 != fold)[0]
        # Withhold complete elemental chains both for evaluation and calibration.
        val_fold=(fold+1)%5
        tr=np.where((z%5!=fold)&(z%5!=val_fold))[0]
        va=np.where(z%5==val_fold)[0];te=np.where(z%5==fold)[0]
        va=np.random.default_rng(5100+fold).permutation(va);tu,ca=np.array_split(va,2)
        yield f'heldout_element_fold{fold}',x[tr],y[tr],x[tu],y[tu],x[ca],y[ca],x[te],y[te]
    tr=np.where(z<=70)[0];va=np.where((z>70)&(z<=82))[0];te=np.where(z>82)[0]
    va=np.random.default_rng(510).permutation(va);tu,ca=np.array_split(va,2)
    yield 'proton_extrapolation',x[tr],y[tr],x[tu],y[tu],x[ca],y[ca],x[te],y[te]
    tr=[];va=[];te=[]
    for v in np.unique(z):
        ids=np.where(z==v)[0];ids=ids[np.argsort(x[ids,1])]
        if len(ids)<8:tr.extend(ids);continue
        tr.extend(ids[:-4]);va.extend(ids[-4:-2]);te.extend(ids[-2:])
    va=np.random.default_rng(511).permutation(va);tu,ca=np.array_split(va,2)
    yield 'neutron_frontier',x[tr],y[tr],x[tu],y[tu],x[ca],y[ca],x[te],y[te]
    for seed in CONFIG['historical_seeds']:
        order=np.random.default_rng(seed).permutation(len(xo))
        tr,tu,ca=np.split(order,[int(.8*len(xo)),int(.9*len(xo))])
        yield f'historical_{seed}',xo[tr],yo[tr],xo[tu],yo[tu],xo[ca],yo[ca],x[later],y[later]

def feature(x,parity=False):
    if not parity:return x
    return np.column_stack([x,CONFIG['parity_coordinate_scale']*(-1.)**x[:,0],
                           CONFIG['parity_coordinate_scale']*(-1.)**x[:,1]])

def kernel(x,z,length,family):
    d2=cdist(x,z,'sqeuclidean')/length**2
    if family=='Matern32':
        q=np.sqrt(3*d2)
        return (1+q)*np.exp(-q)
    return np.exp(-.5*d2)

class GP:
    def __init__(self,x,r,length,noise,family='RBF',parity=False):
        self.x=feature(x,parity);self.length=length;self.family=family;self.parity=parity
        self.signal=max(float(np.std(r)),.1)
        self.nugget=self.signal*noise
        K=self.signal**2*kernel(self.x,self.x,length,family)
        K.flat[::len(K)+1]+=self.nugget**2+1e-10
        self.factor=cho_factor(K,lower=True,check_finite=True)
        self.alpha=cho_solve(self.factor,r,check_finite=True)
    def predict(self,x,uncertainty=False):
        K=self.signal**2*kernel(feature(x,self.parity),self.x,self.length,self.family)
        mean=K@self.alpha
        if not uncertainty:return mean,None
        v=solve_triangular(self.factor[0],K.T,lower=True,check_finite=False)
        variance=self.signal**2-np.sum(v*v,axis=0)+self.nugget**2
        return mean,np.sqrt(np.maximum(variance,1e-10))

def interval_quantile(scores,alpha=.1):
    scores=np.sort(np.asarray(scores))
    rank=math.ceil((len(scores)+1)*(1-alpha))
    if rank>len(scores):return float('inf')
    return float(scores[rank-1])

def fitted_models(x,y,xtune,ytune):
    start=time.perf_counter()
    coef=np.linalg.lstsq(pilot.semf(x),y,rcond=None)[0]
    shellcoef=np.linalg.lstsq(pilot.semf(x,True),y,rcond=None)[0]
    residual=y-pilot.semf(x)@coef
    base=lambda q:pilot.semf(q)@coef
    yield 'SEMF',lambda q:(base(q),np.ones(len(q))),dict(fit_seconds=time.perf_counter()-start)
    yield 'SEMF_known_shells',lambda q:(pilot.semf(q,True)@shellcoef,np.ones(len(q))),dict(fit_seconds=0.)
    start=time.perf_counter();scores=[]
    for length in CONFIG['length_scales']:
        pred=base(xtune)+pilot.residual_predict(x,residual,xtune,length)
        scores.append((pilot.score(ytune,pred)['rmse'],length))
    length=min(scores)[1]
    def smoother(q):
        return base(q)+pilot.residual_predict(x,residual,q,length),np.ones(len(q))
    yield 'Gaussian_residual_smoother',smoother,dict(length=length,fit_seconds=time.perf_counter()-start,
        validation_options=[dict(rmse=s,length=l) for s,l in scores])
    for family,parity,name in [('RBF',False,'GP_RBF'),('Matern32',False,'GP_Matern32'),('RBF',True,'GP_RBF_parity')]:
        start=time.perf_counter();options=[]
        for length in CONFIG['length_scales']:
            for noise in CONFIG['noise_fractions']:
                gp=GP(x,residual,length,noise,family,parity)
                mean,_=gp.predict(xtune)
                risk=pilot.score(ytune,base(xtune)+mean)['rmse']
                options.append((risk,length,noise))
        best=min(options)
        gp=GP(x,residual,best[1],best[2],family,parity)
        def predictor(q,gp=gp):
            r,sd=gp.predict(q,uncertainty=True)
            return base(q)+r,sd
        yield name,predictor,dict(length=best[1],noise_fraction=best[2],fit_seconds=time.perf_counter()-start,
            validation_options=[dict(rmse=a,length=b,noise_fraction=c) for a,b,c in options])

def main():
    out=ROOT/'results/publication';out.mkdir(exist_ok=True)
    (out/'protocol.json').write_text(json.dumps(CONFIG,indent=2))
    (out/'protocol.sha256').write_text(hashlib.sha256((out/'protocol.json').read_bytes()).hexdigest()+'\n')
    rows=[];predictions=[];selections=[];all_splits=[];started=time.perf_counter()
    for tag,x,y,xtu,ytu,xca,yca,xte,yte in partitions():
        ids=lambda a:[f'{int(p)}:{int(n)}' for p,n in a]
        sets=[set(ids(v)) for v in [x,xtu,xca,xte]]
        assert all(not sets[i]&sets[j] for i in range(4) for j in range(i))
        all_splits.append(dict(split=tag,training_ids=ids(x),tuning_ids=ids(xtu),calibration_ids=ids(xca),test_ids=ids(xte)))
        nearest_ca=np.sqrt(cdist(xca,x,'sqeuclidean').min(axis=1))
        nearest_te=np.sqrt(cdist(xte,x,'sqeuclidean').min(axis=1))
        for name,predict,meta in fitted_models(x,y,xtu,ytu):
            t=time.perf_counter();cal,csd=predict(xca);prediction,sd=predict(xte);infer=time.perf_counter()-t
            scale_ca=1+nearest_ca;scale_te=1+nearest_te
            qplain=interval_quantile(abs(yca-cal))
            qdistance=interval_quantile(abs(yca-cal)/scale_ca)
            qgp=interval_quantile(abs(yca-cal)/csd)
            intervals={'constant_residual':np.full(len(xte),qplain),
                'distance_scaled_residual':qdistance*scale_te}
            if name.startswith('GP_'):
                intervals.update(gp_predictive_normal=1.6448536269514722*sd,gp_residual_calibrated=qgp*sd)
            metric=pilot.score(yte,prediction)
            for calibration,width in intervals.items():
                rows.append(dict(split=tag,method=name,interval=calibration,**metric,
                    coverage=float(np.mean(abs(prediction-yte)<=width)),mean_full_width=float(2*np.mean(width)),
                    training_n=len(x),tuning_n=len(xtu),calibration_n=len(xca),fit_seconds=meta['fit_seconds'],
                    calibration_and_test_predict_seconds=infer,mean_test_nearest_training_distance=float(nearest_te.mean())))
                predictions.extend(dict(split=tag,method=name,interval=calibration,Z=int(p[0]),N=int(p[1]),
                    true_binding_MeV=float(t),predicted_binding_MeV=float(v),half_width_MeV=float(w),
                    nearest_training_distance=float(d)) for p,t,v,w,d in zip(xte,yte,prediction,width,nearest_te))
            selections.append(dict(split=tag,method=name,**meta))
            pilot.save_csv(out/'metrics.csv',rows);pilot.save_csv(out/'predictions.csv',predictions)
            (out/'selections.json').write_text(json.dumps(selections,indent=2))
            (out/'splits.json').write_text(json.dumps(all_splits,indent=2))
            print(json.dumps(dict(split=tag,method=name,rmse=metric['rmse'],fit_seconds=meta['fit_seconds'])),flush=True)
    metadata=dict(python=sys.version,numpy=np.__version__,platform=platform.platform(),seconds=time.perf_counter()-started,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        protocol_sha256=hashlib.sha256((out/'protocol.json').read_bytes()).hexdigest(),test_targets_used_for_model_selection=False)
    (out/'environment.json').write_text(json.dumps(metadata,indent=2))
    print(json.dumps(dict(interval_cells=len(rows),prediction_rows=len(predictions),metadata=metadata)))

if __name__=='__main__':main()
