"""Machine-readable scope, residuals, candidates and symmetry convergence."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from full_z_spectral import spectral_operator
from broad_scan import weights


def summarize(root):
    out={'accepted_instabilities':[], 'stability_claim':False, 'runs':[], 'files_sha256':{}}
    for folder in sorted(root.iterdir()):
        if not folder.is_dir() or not (folder/'metadata.json').exists():
            continue
        meta=json.loads((folder/'metadata.json').read_text())
        run={'run':folder.name,'L':meta['L'],'h':meta['h'],
             'stationary_residual':meta['stationary_history'][-1],'sectors':{},'calibration':{}}
        d=np.load(folder/'profile.npz')
        fields=[d[key] for key in ('u','w','s')]
        for m in sorted(map(int,meta['sectors'])):
            reports=[json.loads(p.read_text()) for p in sorted(folder.glob(f'm{m}_*.json'))]
            records=[r for report in reports for r in report['eigenvalues']]
            candidates=[dict(r,strategy=report['strategy'],shift=report['shift'])
                for report in reports for r in report['eigenvalues'] if r['candidate']]
            nonsym=[r for r in candidates if max(r['weighted_symmetry_overlaps'].values(),default=0)<.9]
            unique=[]
            for r in records:
                v=complex(r['real'],r['imag'])
                if not any(abs(v-q)<1e-7 for q in unique): unique.append(v)
            run['sectors'][str(m)]={'targets':len(reports), 'complete_targets':sum(r['status']=='complete' for r in reports),
                'returned_pairs':len(records),'unique_eigenvalues_tolerance_1e-7':len(unique),
                'max_real':max([r['real'] for r in records],default=None),
                'returned_real_range':[min([r['real'] for r in records],default=None),max([r['real'] for r in records],default=None)],
                'returned_imag_range':[min([r['imag'] for r in records],default=None),max([r['imag'] for r in records],default=None)],
                'max_first_order_residual':max([r['first_order_residual'] for r in records],default=None),
                'max_quadratic_scaled_residual':max([r['quadratic_scaled_residual'] for r in records],default=None),
                'max_quadratic_absolute_residual':max([r['quadratic_absolute_residual'] for r in records],default=None),
                'bound':meta['sectors'][str(m)]['bound'],'positive_candidates':candidates,
                'nonsymmetry_candidates':nonsym,
                'shift_targets':[{'strategy':r['strategy'],'shift':r['shift'],'status':r['status']} for r in reports]}
            if m in (0,1) and records:
                _,b,_,modes,_=spectral_operator(fields,meta['n'],meta['h'],meta['c'],meta['nu'],meta['lambda_sigma'],m)
                w=weights(meta['n'],meta['h'],m)
                mask=np.ones((meta['n']-(m!=0),2*meta['n']-1),bool)
                mask[-2:,:]=False;mask[:,:2]=False;mask[:,-2:]=False
                ww=w*np.tile(mask.ravel(),4)
                for name,y in modes.items():
                    best=max(records,key=lambda r:r['weighted_symmetry_overlaps'].get(name,0))
                    by=b@y
                    run['calibration'][name]={'sigma':[best['real'],best['imag']],
                        'overlap':best['weighted_symmetry_overlaps'][name],
                        'boundary_norm_fraction':best['boundary_norm_fraction'],
                        'interior_tangent_residual':float(np.sqrt(np.vdot(by,ww*by).real/np.vdot(y,ww*y).real)),
                        'direct_null_relative_l2':meta['sectors'][str(m)]['neutral_residuals'][name]['relative_l2']}
        out['runs'].append(run)
    for p in sorted(root.rglob('*')):
        if p.is_file() and p.name!='convergence.json':
            out['files_sha256'][p.relative_to(root).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
    out['total_targets']=sum(s['targets'] for r in out['runs'] for s in r['sectors'].values())
    out['total_returned_pairs']=sum(s['returned_pairs'] for r in out['runs'] for s in r['sectors'].values())
    out['nonsymmetry_candidate_count']=sum(len(s['nonsymmetry_candidates']) for r in out['runs'] for s in r['sectors'].values())
    return out


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('root',type=Path)
    args=ap.parse_args()
    report=summarize(args.root)
    (args.root/'convergence.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({k:report[k] for k in ('total_targets','total_returned_pairs','nonsymmetry_candidate_count')}))
