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
        n,h=meta['n'],meta['h']
        carrier=fields[2][:n,1:-1]
        peak=np.unravel_index(np.argmax(carrier),carrier.shape)
        run['background']={'carrier_max':float(carrier.max()),
            'carrier_peak_r_z':[float(peak[0]*h),float((peak[1]-n+1)*h)],
            'carrier_cylindrical_norm_squared_without_2pi':float(np.sum(weights(n,h,0)[:carrier.size]*carrier.ravel()**2)),
            'parity_defects':{key:float(abs(f-sign*f[:,::-1]).max()) for key,f,sign in zip(('u','w','s'),fields,(1,-1,1))}}
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
    out['independent_strategy_agreement']={}
    for m,folder in ((0,'L16_h050_independent'),(2,'L16_h050_independent_m23'),(3,'L16_h050_independent_m23')):
        reference=root/'L16_h050'/f'm{m}_schur_00.json'
        independent=root/folder/f'm{m}_full_00.json'
        if reference.exists() and independent.exists():
            values=lambda path:np.array([complex(r['real'],r['imag']) for r in json.loads(path.read_text())['eigenvalues']])
            a,b=values(reference),values(independent)
            if len(a) and len(b):
                distances=abs(a[:,None]-b[None,:])
                out['independent_strategy_agreement'][str(m)]={
                    'symmetric_max_nearest_eigenvalue_distance':float(max(distances.min(axis=0).max(),distances.min(axis=1).max())),
                    'matching':'Nearest eigenvalue in same sector; no algebraic multiplicity certification',
                    'reference':str(reference.relative_to(root)), 'independent':str(independent.relative_to(root))}
    baseline=root.parent/'full_z_calibration.json'
    out['baseline_calibration_reproduction']=[]
    if baseline.exists():
        old=json.loads(baseline.read_text())
        runs={r['run']:r for r in out['runs']}
        for record in old:
            if record['run'] not in runs:
                continue
            for name,mode in record['neutral_modes'].items():
                new=runs[record['run']]['calibration'].get(name)
                if new:
                    out['baseline_calibration_reproduction'].append({'run':record['run'],'mode':name,
                        'old_abs_sigma':mode['abs_sigma'], 'new_abs_sigma':abs(complex(*new['sigma'])),
                        'absolute_magnitude_difference':abs(abs(complex(*new['sigma']))-mode['abs_sigma']),
                        'phase_direct_null_is_primary_gate':name=='carrier_phase'})
    for p in sorted(root.rglob('*')):
        if p.is_file() and p.name!='convergence.json':
            payload=p.read_bytes()
            if p.suffix in ('.json','.md','.py','.txt'):
                payload=payload.replace(b'\r\n',b'\n')
            out['files_sha256'][p.relative_to(root).as_posix()]=hashlib.sha256(payload).hexdigest()
    out['sha256_scheme']='Text (.json/.md/.py/.txt) normalized CRLF to LF; binary files hashed as raw bytes, so Git line-ending conversion does not invalidate text hashes.'
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
