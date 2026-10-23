"""AME2016 targets only for fitting; evaluate newly non-extrapolated AME2020 masses."""
from pathlib import Path
import csv
import importlib.util
import json
import math
import time
import numpy as np

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('elements',ROOT/'analyze.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
CONFIG=dict(seeds=list(range(7300,7310)),train_fraction=.8,tuning_fraction=.1,
            bandwidths=[1.5,3.,6.],test_definition='Non-extrapolated AME2020 entries absent or mass-estimated in AME2016; A>=16, Z,N>0',
            train_values='AME2016 binding energies only, never updated AME2020 targets')

def load_2016():
    rows={}
    for line in (ROOT/'data/mass16.txt').read_text().splitlines():
        if len(line)<72:continue
        try:N=int(line[4:9]);Z=int(line[9:14]);A=int(line[14:19])
        except ValueError:continue
        assert Z+N==A
        binding=m.number(line[52:63])
        rows[(Z,N)]=dict(Z=Z,N=N,A=A,estimated='#' in line[28:72],binding_per_A_keV=binding)
    return rows

def main():
    start=time.perf_counter();old=load_2016();_,current,_=m.load()
    training=[r for r in old.values() if r['Z']>0 and r['N']>0 and r['A']>=16 and not r['estimated'] and r['binding_per_A_keV'] is not None]
    test=[r for r in current if (r['Z'],r['N']) not in old or old[(r['Z'],r['N'])]['estimated']]
    testids={(r['Z'],r['N']) for r in test};assert not testids&{(r['Z'],r['N']) for r in training}
    x=np.array([[r['Z'],r['N']] for r in training],float);y=np.array([r['binding_per_A_keV']*r['A']/1000 for r in training])
    xt=np.array([[r['Z'],r['N']] for r in test],float);yt=np.array([r['binding_per_A_keV']*r['A']/1000 for r in test])
    rows=[];predictions=[];splits=[]
    for seed in CONFIG['seeds']:
        order=np.random.default_rng(seed).permutation(len(x));a=int(.8*len(x));b=int(.9*len(x));tr=order[:a];tu=order[a:b];ca=order[b:]
        coef=np.linalg.lstsq(m.semf(x[tr]),y[tr],rcond=None)[0]
        shell=np.linalg.lstsq(m.semf(x[tr],True),y[tr],rcond=None)[0]
        residual=y[tr]-m.semf(x[tr])@coef
        options=[]
        for h in CONFIG['bandwidths']:
            v=m.semf(x[tu])@coef+m.residual_predict(x[tr],residual,x[tu],h)
            options.append((m.score(y[tu],v)['rmse'],h))
        h=min(options)[1]
        pt=m.semf(xt)@coef+m.residual_predict(x[tr],residual,xt,h)
        pc=m.semf(x[ca])@coef+m.residual_predict(x[tr],residual,x[ca],h)
        methods={'SEMF':(m.semf(xt)@coef,m.semf(x[ca])@coef),
                 'SEMF_known_shell_features':(m.semf(xt,True)@shell,m.semf(x[ca],True)@shell),
                 'SEFI_geometric_prototype':(pt,pc),'conventional_matched_kernel':(pt.copy(),pc.copy())}
        for name,(pred,calpred) in methods.items():
            rank=min(len(ca),math.ceil((len(ca)+1)*.9));q=float(np.sort(abs(y[ca]-calpred))[rank-1])
            rows.append(dict(seed=seed,method=name,training_n=len(tr),test_n=len(test),bandwidth=h,
                interval_half_width_MeV=q,empirical_90pct_interval_coverage=float(np.mean(abs(pred-yt)<=q)),**m.score(yt,pred)))
            predictions.extend(dict(seed=seed,method=name,Z=r['Z'],N=r['N'],truth_binding_MeV=float(t),prediction_binding_MeV=float(p)) for r,t,p in zip(test,yt,pred))
        splits.append(dict(seed=seed,training_ids=[f'{int(x[i,0])}:{int(x[i,1])}' for i in tr],
            tuning_ids=[f'{int(x[i,0])}:{int(x[i,1])}' for i in tu],calibration_ids=[f'{int(x[i,0])}:{int(x[i,1])}' for i in ca]))
    out=ROOT/'results/historical';out.mkdir(exist_ok=True)
    m.save_csv(out/'metrics.csv',rows);m.save_csv(out/'predictions.csv',predictions)
    (out/'config.json').write_text(json.dumps(CONFIG,indent=2));(out/'splits.json').write_text(json.dumps(splits,indent=2))
    audit=dict(ame2016_rows=len(old),eligible_2016_training_pool=len(training),newly_nonestimated_2020_test=len(test),
        training_test_overlap=0,seconds=time.perf_counter()-start)
    (out/'audit.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit))
    print(json.dumps({name:dict(mean_rmse=float(np.mean([r['rmse'] for r in rows if r['method']==name])),
        mean_coverage=float(np.mean([r['empirical_90pct_interval_coverage'] for r in rows if r['method']==name]))) for name in methods},indent=2))

if __name__=='__main__':main()
