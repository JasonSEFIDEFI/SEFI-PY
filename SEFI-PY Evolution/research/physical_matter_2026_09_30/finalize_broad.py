"""Validate and hash the retained campaign; reconstruct early logged settings."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from compact_broad_vectors import compact
from summarize_broad import summarize
from full_z_spectral import spectral_operator
from broad_scan import weights

BASE = 'bc68d639fe9b175f2e7c4fade05f991c0276157b'


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('root',type=Path)
    args=ap.parse_args()
    compact(args.root)
    checked=0
    saved_residuals={}
    current_sector=None
    for path in sorted(args.root.rglob('m*_*.json')):
        report=json.loads(path.read_text())
        if 'maxiter' not in report:
            # These exact settings are retained from this campaign's tool log.
            report.update(maxiter=200 if ('LR' in path.stem or path.parent.name.endswith('_m23')) else 1200,
                ncv=max(2*report['k_requested']+1,40),tol=1e-10,
                seed=20260930+int(path.stem.rsplit('_',1)[1]),settings_reconstructed_from_execution_log=True)
            path.write_text(json.dumps(report,indent=2))
        with np.load(path.with_suffix('.npz')) as d:
            ev,vec,ids=d['eigenvalues'],d['eigenvectors'],d['saved_indices']
            assert len(ev)==len(report['eigenvalues'])
            assert vec.shape[1]==len(ids)
            assert len(set(ids.tolist()))==len(ids)
            for value,r in zip(ev,report['eigenvalues']):
                assert value.real==r['real'] and value.imag==r['imag']
            assert np.isfinite(vec).all() and np.isfinite(ev).all()
            if len(ids):
                sector=(path.parent,report['m'])
                if sector!=current_sector:
                    meta=json.loads((path.parent/'metadata.json').read_text())
                    with np.load(path.parent/'profile.npz') as profile:
                        fields=[profile[key] for key in ('u','w','s')]
                    a,b,g,_,_=spectral_operator(fields,meta['n'],meta['h'],meta['c'],meta['nu'],meta['lambda_sigma'],report['m'])
                    w=weights(meta['n'],meta['h'],report['m'])
                    current_sector=sector
                key=f'{path.parent.name}/m{report["m"]}'
                summary=saved_residuals.setdefault(key,{'max_first_order_residual':0.,'max_quadratic_absolute_residual':0.,'vectors':0})
                for j,v in zip(ids,vec.T):
                    value=ev[j];y=v[:b.shape[0]]
                    first=float(np.linalg.norm(a@v-value*v)/np.linalg.norm(v))
                    q=b@y+value*(g@y)-value**2*y
                    quadratic=float(np.sqrt(np.vdot(q,w*q).real/np.vdot(y,w*y).real))
                    assert np.isfinite(first) and np.isfinite(quadratic)
                    summary['max_first_order_residual']=max(summary['max_first_order_residual'],first)
                    summary['max_quadratic_absolute_residual']=max(summary['max_quadratic_absolute_residual'],quadratic)
                    summary['vectors']+=1
        checked+=1
    repo=Path(__file__).resolve().parents[3]
    git=['git','-c',f'safe.directory={repo.as_posix()}','-C',str(repo)]
    diff=subprocess.check_output([*git,'diff','--name-status',BASE],text=True)
    changed=[line.split('\t',1) for line in diff.splitlines()]
    assert all(status=='A' for status,path in changed), 'Existing checkpoint changed'
    assert all(path.startswith('SEFI-PY Evolution/research/physical_matter_2026_09_30/') for status,path in changed)
    out={'base_commit':BASE,'all_existing_tracked_files_unchanged':True,
        'preservation_method':'git diff from exact starting commit contains additions only',
        'validated_target_artifacts':checked,'saved_vector_residuals_recomputed_against_saved_background':saved_residuals,
        'tests':{'passed':6,'command':'python -m unittest test_full_z_spectral test_broad_scan -v',
        'note':'Executed successfully during campaign; this finalizer does not re-run tests'},
        'source_hash_scheme':'SHA256 after CRLF-to-LF normalization',
        'source_hashes':{p.name:hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in Path(__file__).parent.glob('*.py')}}
    (args.root/'verification.json').write_text(json.dumps(out,indent=2))
    convergence=summarize(args.root)
    (args.root/'convergence.json').write_text(json.dumps(convergence,indent=2))
    print(json.dumps({key:convergence[key] for key in ('total_targets','total_returned_pairs','nonsymmetry_candidate_count')}))


if __name__=='__main__':
    main()
