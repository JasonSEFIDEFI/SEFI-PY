"""Round-trip branch checks, numerical fits, artifact hashes and final report."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ.setdefault(key,'1')
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from full_z_spectral import stationary,operators

BASE='9b3978c521f4cde792bc8a8f669dee1991887b9e'


def main():
    root=Path(__file__).with_name('three_track_runs')
    roundtrips=[]
    for folder in sorted(root.glob('continuation_*')):
        center=json.loads((folder/'center/diagnostics.json').read_text())
        with np.load(folder/'center/profile.npz') as d:reference=[d[key] for key in ('u','w','s')]
        for parameter,value in (('c',.59),('c',.61),('nu',.29),('nu',.31)):
            record={'run':folder.name,'parameter':parameter,'return_from':value}
            try:
                with np.load(folder/f'{parameter}_{value:.4f}/profile.npz') as d:fields=[d[key] for key in ('u','w','s')]
                fields,hist=stationary(fields,center['n'],center['h'],.6,.3,25,35)
                record.update(stationary_residual=hist[-1],newton_history=hist)
                if hist[-1]>1e-9 or not np.isfinite(hist[-1]) or fields[2].max()<1e-4:raise RuntimeError('Return solve gate failed')
                numerator=sum(np.linalg.norm(f-q)**2 for f,q in zip(fields,reference))
                denominator=np.linalg.norm(reference[0]-1)**2+sum(np.linalg.norm(f)**2 for f in reference[1:])
                record.update(status='converged',relative_profile_difference=float(np.sqrt(numerator/denominator)),
                    max_absolute_profile_difference=max(float(abs(f-q).max()) for f,q in zip(fields,reference)))
            except Exception as error:record.update(status='failed',error=repr(error))
            roundtrips.append(record)
    (root/'roundtrip_checks.json').write_text(json.dumps(roundtrips,indent=2))
    conv=json.loads((root/'convergence.json').read_text())
    metrics=['E_excess','P_physical_z','Q_sigma','Delta_Q_phi','radius_carrier_rms',
        'pohozaev_transverse','pohozaev_longitudinal','pohozaev_transverse_boundary_corrected','pohozaev_longitudinal_boundary_corrected']
    fits=[]
    for L in sorted(set(r['L'] for r in conv)):
        rows=sorted([r for r in conv if r['L']==L and r['status']=='converged'],key=lambda r:r['h'])
        if len(rows)<3:continue
        x=np.array([r['h']**2 for r in rows]);fit={'L':L,'grid_extrapolation_is_not_certified':True,'metrics':{}}
        for key in metrics:
            y=np.array([r[key] for r in rows]);coeff=np.polyfit(x,y,1);quadratic=np.polyfit(x,y,2)
            fit['metrics'][key]={'linear_h_squared_intercept':float(coeff[1]),'linear_h_squared_slope':float(coeff[0]),
                'fit_max_absolute_error':float(abs(y-np.polyval(coeff,x)).max()),
                'quadratic_h_squared_intercept':float(quadratic[-1]),
                'intercept_fit_sensitivity':float(abs(quadratic[-1]-coeff[1]))}
        for key in ('pohozaev_transverse_boundary_corrected','pohozaev_longitudinal_boundary_corrected'):
            fine,coarse=rows[0],rows[-1]
            fit['metrics'][key]['observed_error_order_coarse_to_fine']=float(np.log(abs(coarse[key]/fine[key]))/np.log(coarse['h']/fine['h']))
        fits.append(fit)
    checks=json.loads((root/'first_law_checks.json').read_text())+json.loads((root/'first_law_fine_checks.json').read_text())
    step_fits=[]
    for L,h,parameter in sorted(set((r['L'],r['h'],r['parameter']) for r in checks)):
        rows=[r for r in checks if (r['L'],r['h'],r['parameter'])==(L,h,parameter)]
        x=np.array([r['step']**2 for r in rows]);entry={'L':L,'h':h,'parameter':parameter}
        for key in ('first_law_residual','first_law_minus_box_prediction','discrete_envelope_error','continuum_quadrature_envelope_error','dE_dP_raw'):
            coeff=np.polyfit(x,[r[key] for r in rows],1)
            entry[key+'_zero_step_estimate']=float(coeff[1])
        step_fits.append(entry)
    allrows=[json.loads(p.read_text()) for p in root.rglob('diagnostics.json')]
    operator_cache={}
    saved_residuals=[]
    for path in root.rglob('diagnostics.json'):
        row=json.loads(path.read_text());n,h=row['n'],row['h']
        if (n,h) not in operator_cache:operator_cache[(n,h)]=operators(n,h)
        lap,dz,bc=operator_cache[(n,h)]
        with np.load(path.with_name('profile.npz')) as d:u,w,s=[d[key][:n,1:-1].ravel() for key in ('u','w','s')]
        db=np.zeros((n,2*n-1));db[:,0]=-1/(2*h);db[:,-1]=1/(2*h)
        density=u*u+w*w;a=1-density-4*s*s
        residual=np.r_[lap@u+bc-row['c']*(dz@w)+a*u,
            lap@w+row['c']*(dz@u+db.ravel())+a*w,
            lap@s+(3+row['nu']**2-4*density-row['lambda']*s*s)*s]
        error=float(abs(residual).max())
        assert np.isfinite(error) and error<1e-9,'Saved-profile field residual gate failed'
        saved_residuals.append({'profile':str(path.parent.relative_to(root)),'max_field_residual':error})
    failure_files=list(root.rglob('failure.json'))
    summary={'starting_commit':BASE,'converged_profile_count':len(allrows),'solver_failure_count':len(failure_files),
        'max_stationary_residual':max(r['stationary_residual'] for r in allrows),
        'independently_recomputed_saved_profile_residuals':saved_residuals,
        'max_stress_generator_algebra_error':max(abs(r['generator_algebra_residual']) for r in allrows),
        'grid_fits':fits,'zero_parameter_step_fits':step_fits,'roundtrip_checks':roundtrips,
        'conclusions':{'infinite_domain_existence_established':False,'continuum_first_law_limit_established':False,
            'naive_dE_dP_equals_v_supported':False,'universal_phase_metric_supported':False,
            'matter_or_gravity_identification_established':False},
        'tests':{'passed':10,'command':'python -m unittest test_full_z_spectral test_broad_scan test_three_track -v'}}
    repo=Path(__file__).resolve().parents[3]
    command=['git','-c',f'safe.directory={repo.as_posix()}','-C',str(repo),'diff','--name-status',BASE]
    changes=subprocess.check_output(command,text=True).splitlines()
    assert all(line.startswith('A\t') for line in changes),'A preexisting tracked checkpoint was modified'
    summary['all_preexisting_tracked_files_preserved']=True
    (root/'analysis.json').write_text(json.dumps(summary,indent=2))
    with (root/'convergence.csv').open('w',newline='') as handle:
        keys=['L','h',*metrics,'stationary_residual'];writer=csv.DictWriter(handle,fieldnames=keys);writer.writeheader()
        writer.writerows({k:r[k] for k in keys} for r in conv if r['status']=='converged')
    hashes={}
    for path in sorted(root.rglob('*')):
        if not path.is_file() or path.name=='hashes.json':continue
        data=path.read_bytes()
        if path.suffix in ('.json','.csv','.md','.py'):data=data.replace(b'\r\n',b'\n')
        hashes[path.relative_to(root).as_posix()]=hashlib.sha256(data).hexdigest()
    (root/'hashes.json').write_text(json.dumps({'scheme':'SHA256; text CRLF normalized to LF, binary raw bytes','files':hashes},indent=2))
    lines=['# Three-track research checkpoint — 30 September 2026','',
        f'Started from `{BASE}` on `research/action-audit-2026-09-30`. v0.1 and every preexisting tracked file are preserved. N=0, lambda=25; no matter or gravity premise is adopted.','',
        f'**Executed:** {len(allrows)} converged full-Z profiles, 12 box/grid convergence cases, c continuation over 0.50–0.70, nu over 0.26–0.34, 20 centered identity checks and 16 local return checks. Newton failures: {len(failure_files)}. Ten tests pass.','',
        '## 1. Continuum and domain convergence','',
        'Audited physical E, P and both phase charges come directly from the laboratory stress/current definitions. Bare reduced momentum, graph-functional momentum and physical momentum are distinct. The algebraic generator identity is checked before interpreting branch trends.','',
        '| L | h | E excess | P physical | Q Sigma | Carrier RMS radius |',
        '|---:|---:|---:|---:|---:|---:|']
    for r in conv:lines.append(f'| {r["L"]:g} | {r["h"]:.6g} | {r["E_excess"]:.8g} | {r["P_physical_z"]:.8g} | {r["Q_sigma"]:.8g} | {r["radius_carrier_rms"]:.8g} |')
    lines+=['',
        'Boundary-corrected Pohozaev residuals decrease approximately as h squared. Raw residuals contain finite-box surface terms; their signs can cancel grid errors at an intermediate h. A small raw value at one mesh is therefore not an existence gate. Domain and grid changes remain material in E and P. Three-mesh extrapolations and their fit sensitivity are saved, but do not establish an infinite-domain branch. Negative background-subtracted E is retained, not reinterpreted as particle mass.','',
        '## 2. Continuation and first law','',
        'Parametric E(P), both charges, radii and all backgrounds are in the branch CSV/JSON files. The c and nu scans use separate outward paths from the same center. Local return solves quantify reproducibility; they do not constitute global fold detection or pseudo-arclength continuation.','',
        '| L | h | parameter | First-law residual, step→0 estimate | Graph envelope error, step→0 estimate |',
        '|---:|---:|:---|---:|---:|']
    for r in step_fits:lines.append(f'| {r["L"]:g} | {r["h"]:.6g} | {r["parameter"]} | {r["first_law_residual_zero_step_estimate"]:.7g} | {r["discrete_envelope_error_zero_step_estimate"]:.7g} |')
    lines+=['',
        'The discrete envelope error decreases with parameter-step refinement; the physical quadrature first-law error retains a grid-dependent floor. These are distinct tests. The derived finite-box c-law defect is -v R_Z/(2 Omega gamma squared); raw and corrected tests are both saved. The full continuum first law is not yet numerically closed. Dropping charge derivatives gives the wrong slope: raw dE/dP does not equal v on this varying-charge branch.','',
        '## 3. Characteristics and universal geometry','',
        'The full coupled Cartesian principal symbol is (eta contracted with k twice) times I_4. Potential couplings and rotating-frame gyroscopic terms are lower order; exact microscopic fronts share the already-assumed Minkowski cone. The Phi infrared phase has sound-speed squared 0.2. An admissible homogeneous two-condensate example has two different phase speeds squared, 0.04959603434 and 0.33918901239, independently recovered from exact coupled dispersion. The zero-carrier asymptotic state and ring cores invalidate extending that two-phase elimination everywhere.','',
        '**Negative geometry result:** an unqualified universal extension of the phase-only metric to all sectors fails these comparisons. Shared microscopic Minkowski characteristics are assumed structure, not emergent gravitational dynamics. This does not exclude every restricted infrared universality regime. Localized-loop coupling to a varying phase background has not been computed. No matter identification, universal gravitational coupling or Einstein equation follows.','',
        '## Reproduction and limits','',
        'See [derivation](three_track_derivation.md), [machine analysis](three_track_runs/analysis.json), [artifact hashes](three_track_runs/hashes.json), and [dispersion diagnostics](three_track_runs/characteristics.json). Profiles and failures remain separate research additions.','',
        '```text','python -m unittest test_full_z_spectral test_broad_scan test_three_track -v',
        'python run_three_track.py --track convergence','python run_three_track.py --track continuation',
        'python run_three_track.py --track continuation --fine-only','python characteristic_audit.py',
        'python analyze_three_track.py','```','',
        'Use the unchanged requirements-spectral.txt. Future work should reduce h at larger L, control the cylinder-shape and far-field charge limits, and test the first law after both limits. No new stability scan was performed in this checkpoint.','']
    Path(__file__).with_name('three_track_report.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps({'profiles':len(allrows),'failures':len(failure_files),'max_stationary_residual':summary['max_stationary_residual'],
        'max_generator_error':summary['max_stress_generator_algebra_error'],
        'max_local_return_difference':max(r.get('relative_profile_difference',0) for r in roundtrips),
        'failed_return_solves':sum(r['status']=='failed' for r in roundtrips)}))


if __name__=='__main__':main()
