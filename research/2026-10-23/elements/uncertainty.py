from pathlib import Path
import csv
import json
import numpy as np
ROOT=Path(__file__).resolve().parent

def main():
    rows=list(csv.DictReader((ROOT/'results/mass_predictions.csv').open()))
    historical=ROOT/'results/historical/predictions.csv'
    if historical.exists():
        for r in csv.DictReader(historical.open()):
            r['split']='historical_'+r['seed'];rows.append(r)
    methods=sorted({r['method'] for r in rows});splits=sorted({r['split'] for r in rows})
    results=[]
    for tag in splits:
        baseline=[r for r in rows if r['split']==tag and r['method']=='SEMF']
        z=np.array([int(r['Z']) for r in baseline]);groups=np.unique(z)
        base=np.array([(float(r['prediction_binding_MeV'])-float(r['truth_binding_MeV']))**2 for r in baseline])
        base_sum=np.array([base[z==g].sum() for g in groups]);counts=np.array([(z==g).sum() for g in groups])
        ix=np.random.default_rng(824).integers(len(groups),size=(2000,len(groups)))
        for name in methods:
            selected=[r for r in rows if r['split']==tag and r['method']==name]
            assert [(r['Z'],r['N']) for r in selected]==[(r['Z'],r['N']) for r in baseline]
            error=np.array([(float(r['prediction_binding_MeV'])-float(r['truth_binding_MeV']))**2 for r in selected])
            sums=np.array([error[z==g].sum() for g in groups])
            delta=np.sqrt(sums[ix].sum(axis=1)/counts[ix].sum(axis=1))-np.sqrt(base_sum[ix].sum(axis=1)/counts[ix].sum(axis=1))
            results.append(dict(split=tag,method=name,proton_clusters=len(groups),
                rmse_delta_vs_semf=float(np.sqrt(error.mean())-np.sqrt(base.mean())),
                cluster_bootstrap_ci95=np.quantile(delta,[.025,.975]).tolist()))
    (ROOT/'results/uncertainty.json').write_text(json.dumps(results,indent=2));print('Saved',len(results),'cluster bootstrap comparisons')

if __name__=='__main__':main()
