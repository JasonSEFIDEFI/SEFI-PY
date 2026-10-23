"""Exploratory, leakage-controlled isotope prediction from public evaluated data.

No outputs are claims of a new element, nuclear existence, or SEFI dynamics.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import platform
import re
import sys
import time
import numpy as np
from scipy.spatial.distance import cdist

ROOT=Path(__file__).resolve().parent
CONFIG=dict(seed_range=list(range(4200,4210)),minimum_A=16,
    kernel_bandwidths=[1.5,3.,6.],kernel_prior_weight=1.,shell_width=3.,
    mass_primary_endpoint='Total binding-energy RMSE in MeV on held-out nuclei',
    lifetime_endpoint='MAE of log10(half-life / 1 second), finite measured radioactive ground states only',
    proton_magic=[2,8,20,28,50,82],neutron_magic=[2,8,20,28,50,82,126],
    status='Exploratory reconstruction of a lost original map; conventional models and matched geometric re-expression.')

def number(field):
    value=field.strip()
    if not value or '*' in value:return None
    return float(value.replace('#','.'))

def load():
    ame={}; nubase={};rejected=[]
    for i,line in enumerate((ROOT/'data/mass_1.mas20.txt').read_text().splitlines(),1):
        if len(line)<78:continue
        try:N=int(line[4:9]); Z=int(line[9:14]); A=int(line[14:19])
        except ValueError:continue
        if A!=Z+N:raise ValueError(('bad AME counts',i))
        b=number(line[54:67]);u=number(line[68:78])
        row=dict(Z=Z,N=N,A=A,symbol=line[20:23].strip(),binding_per_A_keV=b,
                 binding_unc_per_A_keV=u,mass_excess_keV=number(line[28:42]),
                 mass_unc_keV=number(line[42:54]),mass_estimated='#' in line[28:78],ame_line=i)
        key=(Z,N)
        if key in ame:raise ValueError(('duplicate AME',key))
        ame[key]=row
    units={'ys':1e-24,'zs':1e-21,'as':1e-18,'fs':1e-15,'ps':1e-12,'ns':1e-9,'us':1e-6,'ms':1e-3,
           's':1,'m':60,'h':3600,'d':86400,'y':365.25*86400,'ky':1e3*365.25*86400,
           'My':1e6*365.25*86400,'Gy':1e9*365.25*86400,'Ty':1e12*365.25*86400,'Py':1e15*365.25*86400,
           'Ey':1e18*365.25*86400,'Zy':1e21*365.25*86400,'Yy':1e24*365.25*86400}
    for i,line in enumerate((ROOT/'data/nubase_4.mas20.txt').read_text().splitlines(),1):
        if line.startswith('#') or len(line)<80:continue
        A=int(line[0:3]); Z=int(line[4:7]);state=line[7:8]
        if state!='0':continue
        field=line[69:78].strip();unit=line[78:80].strip()
        finite=None;status='missing'
        if field=='stbl':status='no_observed_decay'
        elif field=='p-unst':status='particle_unstable'
        elif '#' in field:status='estimated'
        elif any(c in field for c in '<>~'):status='bounded_or_approximate'
        elif re.fullmatch(r'[0-9.+\-Ee]+',field or ' '):
            if unit not in units:raise ValueError(('unknown lifetime unit',i,unit))
            finite=float(field)*units[unit];status='finite_evaluated'
            if finite<=0:raise ValueError(('nonpositive lifetime',i))
        try:year=int(line[114:118])
        except ValueError:year=None
        nubase[(Z,A-Z)]=dict(half_life_seconds=finite,lifetime_status=status,
            lifetime_raw=field,lifetime_unit=unit,spin_parity=line[88:102].strip(),
            discovery_year=year,decay_modes=line[119:209].strip(),nubase_line=i)
    rows=[]
    for key,row in sorted(ame.items()):
        row.update(nubase.get(key,{}));rows.append(row)
    audit=dict(ame_rows=len(ame),nubase_ground_rows=len(nubase),joined_rows=len(rows),
        joined_nubase=sum('lifetime_status' in row for row in rows),
        estimated_mass_rows=sum(row['mass_estimated'] for row in rows),
        min_Z=min(row['Z'] for row in rows),max_Z=max(row['Z'] for row in rows),
        lifetime_status_counts={s:sum(row.get('lifetime_status')==s for row in rows) for s in sorted({r.get('lifetime_status','unmatched') for r in rows})})
    selected=[r for r in rows if r['Z']>0 and r['N']>0 and r['A']>=CONFIG['minimum_A'] and
              not r['mass_estimated'] and r['binding_per_A_keV'] is not None]
    audit['mass_analysis_rows']=len(selected)
    return rows,selected,audit

def semf(x, shells=False):
    Z=x[:,0];N=x[:,1];A=Z+N
    pairing=(((-1.)**Z+(-1.)**N)/2)/np.sqrt(A)
    terms=[A,-A**(2/3),-Z*(Z-1)/A**(1/3),-(N-Z)**2/A,pairing]
    if shells:
        dz=np.min(abs(Z[:,None]-np.array(CONFIG['proton_magic'])),axis=1)
        dn=np.min(abs(N[:,None]-np.array(CONFIG['neutron_magic'])),axis=1)
        p=np.exp(-.5*(dz/CONFIG['shell_width'])**2);n=np.exp(-.5*(dn/CONFIG['shell_width'])**2)
        terms += [p,n,p*n]
    return np.column_stack(terms)

def residual_predict(xtrain, residual, x, bandwidth):
    w=np.exp(-cdist(x,xtrain,'sqeuclidean')/(2*bandwidth**2))
    return (w@residual)/(w.sum(axis=1)+CONFIG['kernel_prior_weight'])

def split_masks(rows):
    n=len(rows);Z=np.array([r['Z'] for r in rows]);N=np.array([r['N'] for r in rows]);years=np.array([r.get('discovery_year') or 9999 for r in rows])
    for seed in CONFIG['seed_range']:
        order=np.random.default_rng(seed).permutation(n);parts=np.array_split(order,[int(.6*n),int(.8*n)])
        yield f'random_{seed}',*parts
    yield 'proton_extrapolation',np.where(Z<=70)[0],np.where((Z>70)&(Z<=82))[0],np.where(Z>82)[0]
    train=[];valid=[];test=[]
    for z in np.unique(Z):
        ids=np.where(Z==z)[0];ids=ids[np.argsort(N[ids])]
        if len(ids)<8:train.extend(ids);continue
        train.extend(ids[:-4]);valid.extend(ids[-4:-2]);test.extend(ids[-2:])
    yield 'neutron_frontier',np.array(train),np.array(valid),np.array(test)
    yield 'discovery_year_proxy',np.where(years<=2000)[0],np.where((years>2000)&(years<=2010))[0],np.where((years>2010)&(years<=2020))[0]

def score(truth,pred):
    delta=pred-truth
    return dict(n=len(truth),rmse=float(np.sqrt(np.mean(delta**2))),mae=float(np.mean(abs(delta))),bias=float(delta.mean()))

def half_life(rows, tr,va,te,tag):
    x=np.array([[r['Z'],r['N']] for r in rows],dtype=float)
    finite=np.array([r.get('half_life_seconds') is not None for r in rows])
    train=tr[finite[tr]];valid=va[finite[va]];test=te[finite[te]]
    stable=np.array([i for i in tr if rows[i].get('lifetime_status')=='no_observed_decay'])
    if min(len(train),len(valid),len(test),len(stable))<5:return [],[]
    y=np.array([math.log10(r['half_life_seconds']) if r.get('half_life_seconds') else np.nan for r in rows])
    distance=np.sqrt(cdist(x,x[stable],'sqeuclidean').min(axis=1))
    beta=np.linalg.lstsq(np.column_stack((np.ones(len(train)),distance[train])),y[train],rcond=None)[0]
    # Enforce the proposed nonincreasing lifetime-distance law.
    if beta[1]>0:beta=np.array([y[train].mean(),0.])
    def poly(a):
        z=a[:,0]/100;n=a[:,1]/150
        return np.column_stack((np.ones(len(a)),z,n,z*z,z*n,n*n,(-1.)**a[:,0],(-1.)**a[:,1]))
    coeff=np.linalg.lstsq(poly(x[train]),y[train],rcond=None)[0]
    center=y[train].mean()
    candidates=[]
    for h in CONFIG['kernel_bandwidths']:
        v=center+residual_predict(x[train],y[train]-center,x[valid],h)
        candidates.append((score(y[valid],v)['mae'],h))
    h=min(candidates)[1]
    predictions={'training_median':np.full(len(test),np.median(y[train])),
        'distance_stability_hypothesis':beta[0]+beta[1]*distance[test],
        'conventional_polynomial':poly(x[test])@coeff,
        'conventional_kernel':center+residual_predict(x[train],y[train]-center,x[test],h)}
    results=[];raw=[]
    for method,pred in predictions.items():
        results.append(dict(split=tag,method=method,**score(y[test],pred),training_n=len(train),stable_anchors=len(stable),bandwidth=h,fit_distance_slope=float(beta[1])))
        raw.extend(dict(split=tag,method=method,Z=rows[i]['Z'],N=rows[i]['N'],truth_log10_seconds=float(y[i]),prediction_log10_seconds=float(v)) for i,v in zip(test,pred))
    return results,raw

def save_csv(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def main():
    out=ROOT/'results';out.mkdir(exist_ok=True)
    (out/'config.json').write_text(json.dumps(CONFIG,indent=2));start=time.perf_counter()
    allrows,rows,audit=load();save_csv(out/'evaluated_ground_states.csv',allrows)
    (out/'data_audit.json').write_text(json.dumps(audit,indent=2))
    x=np.array([[r['Z'],r['N']] for r in rows],float)
    y=np.array([r['binding_per_A_keV']*r['A']/1000 for r in rows])
    results=[];predictions=[];life_results=[];life_raw=[];splits=[]
    for tag,tr,va,te in split_masks(rows):
        assert not set(tr)&set(va) and not set(tr)&set(te) and not set(va)&set(te)
        if min(len(tr),len(va),len(te))<10:continue
        coef=np.linalg.lstsq(semf(x[tr]),y[tr],rcond=None)[0]
        shellcoef=np.linalg.lstsq(semf(x[tr],True),y[tr],rcond=None)[0]
        residual=y[tr]-semf(x[tr])@coef
        tune=va[:len(va)//2];cal=va[len(va)//2:]
        options=[]
        for h in CONFIG['kernel_bandwidths']:
            v=semf(x[tune])@coef+residual_predict(x[tr],residual,x[tune],h)
            options.append((score(y[tune],v)['rmse'],h))
        bandwidth=min(options)[1]
        smooth=residual_predict(x[tr],residual,x[te],bandwidth)
        # Independent Euclidean implementation of the same Gaussian weights.
        d2=((x[te,None,:]-x[None,tr,:])**2).sum(axis=2)
        w=np.exp(-d2/(2*bandwidth**2));ordinary=(w@residual)/(w.sum(axis=1)+1)
        assert np.allclose(smooth,ordinary,atol=1e-12)
        preds={'SEMF':semf(x[te])@coef,'SEMF_known_shell_features':semf(x[te],True)@shellcoef,
            'SEFI_geometric_prototype':semf(x[te])@coef+smooth,
            'conventional_matched_kernel':semf(x[te])@coef+ordinary}
        calbase=semf(x[cal])@coef
        calsmooth=calbase+residual_predict(x[tr],residual,x[cal],bandwidth)
        calpred={'SEMF':calbase,'SEMF_known_shell_features':semf(x[cal],True)@shellcoef,
                 'SEFI_geometric_prototype':calsmooth,'conventional_matched_kernel':calsmooth}
        for method,pred in preds.items():
            rank=min(len(cal),math.ceil((len(cal)+1)*.9))
            q=float(np.sort(abs(y[cal]-calpred[method]))[rank-1])
            results.append(dict(split=tag,method=method,training_n=len(tr),validation_n=len(va),test_n=len(te),bandwidth=bandwidth,
                calibration_n=len(cal),interval_half_width_MeV=q,empirical_90pct_interval_coverage=float(np.mean(abs(pred-y[te])<=q)),**score(y[te],pred)))
            predictions.extend(dict(split=tag,method=method,Z=rows[i]['Z'],N=rows[i]['N'],A=rows[i]['A'],truth_binding_MeV=float(y[i]),prediction_binding_MeV=float(v)) for i,v in zip(te,pred))
        splits.append(dict(split=tag,training_ids=[f'{rows[i]["Z"]}:{rows[i]["N"]}' for i in tr],
            validation_ids=[f'{rows[i]["Z"]}:{rows[i]["N"]}' for i in va],test_ids=[f'{rows[i]["Z"]}:{rows[i]["N"]}' for i in te],
            tuning_ids=[f'{rows[i]["Z"]}:{rows[i]["N"]}' for i in tune],calibration_ids=[f'{rows[i]["Z"]}:{rows[i]["N"]}' for i in cal],coefficients=coef.tolist(),bandwidth=bandwidth))
        a,b=half_life(rows,tr,va,te,tag);life_results.extend(a);life_raw.extend(b)
        print(tag,{method:score(y[te],p)['rmse'] for method,p in preds.items()},flush=True)
    save_csv(out/'mass_metrics.csv',results);save_csv(out/'mass_predictions.csv',predictions)
    save_csv(out/'lifetime_metrics.csv',life_results);save_csv(out/'lifetime_predictions.csv',life_raw)
    (out/'splits.json').write_text(json.dumps(splits,indent=2))
    metadata=dict(python=sys.version,numpy=np.__version__,platform=platform.platform(),seconds=time.perf_counter()-start,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),new_elements_predicted=0,
        evaluated_isotope_gaps_are_not_existence_predictions=True)
    (out/'environment.json').write_text(json.dumps(metadata,indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    pts=axes[0].scatter(x[:,1],x[:,0],c=y/(x.sum(axis=1)),s=7,cmap='viridis')
    axes[0].set(xlabel='Neutron number N',ylabel='Proton number Z',title='Measured AME2020 ground-state mass subset')
    fig.colorbar(pts,ax=axes[0],label='Binding energy per nucleon (MeV)')
    labels=['random_4200','neutron_frontier','proton_extrapolation'];xx=np.arange(3)
    for j,method in enumerate(['SEMF','SEMF_known_shell_features','SEFI_geometric_prototype','conventional_matched_kernel']):
        vals=[next(r['rmse'] for r in results if r['split']==s and r['method']==method) for s in labels]
        axes[1].bar(xx+(j-1.5)*.2,vals,.2,label=method)
    axes[1].set(xticks=xx,xticklabels=['Random','Neutron frontier','Proton extrapolation'],ylabel='Total binding-energy RMSE (MeV)',title='Held-out errors; prototype and matched baseline coincide')
    axes[1].legend(fontsize=6);fig.tight_layout();fig.savefig(out/'isotope_mapping_and_prediction.png',dpi=160);plt.close(fig)
    print(json.dumps(dict(audit=audit,mass_cells=len(results),lifetime_cells=len(life_results),metadata=metadata),indent=2))

if __name__=='__main__':main()
