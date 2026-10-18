"""Frozen simulation, controls, paired statistics, timing, and physical checks.

Run from this directory: python run_benchmark.py
All outputs are derived from protocol.json; no selection of favorable seeds.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
import csv
import hashlib
import json
import platform
import sys
import time
from pathlib import Path
import numpy as np
from bosonic import metrics, code_overlap, displacement, erasure_kraus, choi
from decoders import covariance_from_calibration, make_decoders
from statistics import interval

ROOT = Path(__file__).resolve().parent


def save_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def noise_covariance(sigma, ratio, angle):
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    return rotation@np.diag([sigma*sigma, (sigma*ratio)**2])@rotation.T


def physical_validation(out):
    angles = np.linspace(0,2*np.pi,9)[:-1]
    betas = np.array([r*np.exp(1j*a) for r in np.linspace(0,3,13) for a in angles])
    rows=[]
    for cutoff in (8,16,32,48):
        max_error=0.; max_unitary=0.
        for beta in betas:
            d=displacement(complex(beta),cutoff)
            max_error=max(max_error,float(np.max(np.abs(d[:2,:2]-code_overlap(beta)))))
            max_unitary=max(max_unitary,float(np.max(np.abs(d.conj().T@d-np.eye(cutoff)))))
        rows.append({"cutoff":cutoff,"max_code_block_error":max_error,
                     "max_unitarity_error":max_unitary,"residual_radius_max":3.0,"points":len(betas)})
    kraus_defect=0.; min_choi=1.
    for beta in betas:
        ks=erasure_kraus(beta)
        kraus_defect=max(kraus_defect,float(np.max(np.abs(sum(k.conj().T@k for k in ks)-np.eye(2)))))
        min_choi=min(min_choi,float(np.linalg.eigvalsh(choi(ks)).min()))
    save_csv(out/'cutoff_convergence.csv',rows)
    physical={"kraus_completeness_max_error":kraus_defect,"minimum_choi_eigenvalue":min_choi,
              "cutoff_48_max_overlap_error":rows[-1]['max_code_block_error'],
              "scores_use_exact_infinite_oscillator_overlaps":True,
              "cutoff_grid_covers_radius_0_to_3":True}
    assert kraus_defect<1e-12 and min_choi>-1e-12 and rows[-1]['max_code_block_error']<1e-10
    return physical


def run(protocol, out):
    out.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter(); n=protocol['trials_per_seed_per_noise']
    sigmas=protocol['noise_sigmas']; seeds=protocol['seeds']
    r=np.eye(2)/(2*protocol['pilot_modes_per_trial'])
    transform=np.asarray(protocol['coordinate_matrix'])
    observations=np.empty((seeds,len(sigmas),n,2)); truth=np.empty_like(observations)
    calibrated=np.empty((seeds,len(sigmas),2,2))
    seed_rows=[]; stress_rows=[]; audit_max={"affine_estimate_error":0.,"nonlinear_estimate_error":0.,"conventional_estimate_error":0.}
    method_order=None; setup_times={}; timing_decoders=None; timing_y=None
    for j,sigma in enumerate(sigmas):
        population_cov=noise_covariance(sigma,protocol['noise_axis_ratio'],protocol['noise_angle_radians'])
        for k in range(seeds):
            # Independent RNG streams: calibration, held-out displacement, measurement, controls.
            streams=np.random.SeedSequence([protocol['first_seed']+k,j]).spawn(4)
            train_rng,noise_rng,measure_rng,control_rng=[np.random.default_rng(s) for s in streams]
            calibration=train_rng.multivariate_normal([0,0],population_cov,protocol['calibration_per_seed_per_noise'])
            fitted=covariance_from_calibration(calibration); calibrated[k,j]=fitted
            beta=noise_rng.multivariate_normal([0,0],population_cov,n)
            # Sufficient statistic of M independent heterodyne probes with known alpha.
            # Each probe has real/imaginary shot-noise variance 1/2; sample mean has R.
            y=beta+measure_rng.multivariate_normal([0,0],r,n)
            truth[k,j]=beta; observations[k,j]=y
            t=time.perf_counter(); decoders=make_decoders(fitted,r,transform,control_rng,n)
            setup_times.setdefault(str(sigma),[]).append(time.perf_counter()-t)
            if method_order is None:method_order=list(decoders)
            predicted={name:decoder.predict(y) for name,decoder in decoders.items()}
            audit_max['affine_estimate_error']=max(audit_max['affine_estimate_error'],float(np.max(np.abs(predicted['Geometry']-predicted['Affine coordinates']))))
            audit_max['nonlinear_estimate_error']=max(audit_max['nonlinear_estimate_error'],float(np.max(np.abs(predicted['Geometry']-predicted['Nonlinear coordinates']))))
            audit_max['conventional_estimate_error']=max(audit_max['conventional_estimate_error'],float(np.max(np.abs(predicted['Geometry']-predicted['Conventional posterior']))))
            for name,estimate in predicted.items():
                scores=metrics(beta-estimate,protocol['accuracy_residual_radius'])
                row={'sigma':sigma,'seed':protocol['first_seed']+k,'method':name,'n_trials':n}
                row.update({metric:float(values.mean()) for metric,values in scores.items()})
                seed_rows.append(row)
            if sigma==protocol['primary_noise_sigma'] and k==0:
                timing_decoders=decoders; timing_y=y
            if sigma==protocol['stress_sigma']:
                shift=np.asarray(protocol['stress_mean_shift'])
                shifted_beta=beta+shift; shifted_y=y+shift
                for name in ('No recovery','ML','Geometry','Conventional posterior'):
                    scores=metrics(shifted_beta-decoders[name].predict(shifted_y),protocol['accuracy_residual_radius'])
                    row={'sigma':sigma,'seed':protocol['first_seed']+k,'method':name,'n_trials':n,'scenario':'unmodeled mean shift'}
                    row.update({metric:float(values.mean()) for metric,values in scores.items()})
                    stress_rows.append(row)
        print(f'Completed sigma={sigma}: {seeds*n:,} held-out trials, all controls.',flush=True)
    save_csv(out/'seed_metrics.csv',seed_rows); save_csv(out/'stress_seed_metrics.csv',stress_rows)
    np.savez_compressed(out/'raw_trials.npz',truth_beta=truth,observations_y=observations,
                        calibration_covariances=calibrated,sigmas=np.asarray(sigmas),measurement_covariance=r)
    aggregate=[]
    for sigma in sigmas:
        ml=[row for row in seed_rows if row['sigma']==sigma and row['method']=='ML']
        geo=[row for row in seed_rows if row['sigma']==sigma and row['method']=='Geometry']
        for name in method_order:
            selected=[row for row in seed_rows if row['sigma']==sigma and row['method']==name]
            row={'sigma':sigma,'method':name,'seeds':seeds,'trials':seeds*n}
            for metric in ('entanglement_fidelity','state_fidelity','survival','accuracy','residual_mse'):
                mean,lo,hi=interval([s[metric] for s in selected],protocol['bootstrap_replicates'],protocol['bootstrap_seed'])
                row.update({metric:mean,metric+'_ci_low':lo,metric+'_ci_high':hi})
            differences=np.array([s['entanglement_fidelity']-m['entanglement_fidelity'] for s,m in zip(selected,ml)])
            effect,lo,hi=interval(differences,protocol['bootstrap_replicates'],protocol['bootstrap_seed'])
            row.update({'paired_delta_fe_vs_ml':effect,'paired_delta_ci_low':lo,'paired_delta_ci_high':hi,
                        'relative_infidelity_reduction_vs_ml':effect/(1-np.mean([m['entanglement_fidelity'] for m in ml]))})
            geffect,glo,ghi=interval([s['entanglement_fidelity']-g['entanglement_fidelity'] for s,g in zip(selected,geo)],protocol['bootstrap_replicates'],protocol['bootstrap_seed'])
            row.update({'paired_delta_fe_vs_geometry':geffect,'paired_delta_geometry_ci_low':glo,'paired_delta_geometry_ci_high':ghi})
            aggregate.append(row)
    save_csv(out/'table_1_all_results.csv',aggregate)
    stress=[]
    stress_ml=[x for x in stress_rows if x['method']=='ML']
    for name in ('No recovery','ML','Geometry','Conventional posterior'):
        selected=[x for x in stress_rows if x['method']==name]
        mean,lo,hi=interval([x['entanglement_fidelity'] for x in selected],protocol['bootstrap_replicates'],protocol['bootstrap_seed'])
        effect,elo,ehi=interval([x['entanglement_fidelity']-m['entanglement_fidelity'] for x,m in zip(selected,stress_ml)],protocol['bootstrap_replicates'],protocol['bootstrap_seed'])
        stress.append({'method':name,'fe':mean,'ci_low':lo,'ci_high':hi,'paired_delta_vs_ml':effect,'delta_ci_low':elo,'delta_ci_high':ehi})
    save_csv(out/'distribution_shift.csv',stress)
    # Decode-only costs: common preparation/scoring excluded, same available y.
    # Median of repeated whole batches. Rotating method order reduces ordering bias.
    times={name:[] for name in method_order}; timing_rng=np.random.default_rng(31021)
    for decoder in timing_decoders.values():decoder.predict(timing_y)
    for rep in range(protocol['runtime_repeats']):
        for name in timing_rng.permutation(method_order):
            t=time.perf_counter_ns(); timing_decoders[name].predict(timing_y)
            times[name].append((time.perf_counter_ns()-t)/1000/len(timing_y))
    runtimes=[]
    for name,values in times.items():
        lo,hi=np.quantile(values,[.025,.975])
        runtimes.append({'method':name,'median_microseconds_per_trial':float(np.median(values)),
                         'timing_p025':float(lo),'timing_p975':float(hi),'repeats':protocol['runtime_repeats'],
                         'batch':len(timing_y),'timing_range_is_not_statistical_ci':True})
    save_csv(out/'runtime.csv',runtimes)
    timing_lookup={row['method']:row for row in runtimes}
    for row in aggregate:
        timing=timing_lookup[row['method']]
        row.update({'decode_median_us_per_trial_at_sigma_0_2':timing['median_microseconds_per_trial'],
                    'runtime_p025':timing['timing_p025'],'runtime_p975':timing['timing_p975'],
                    'cost_scope':'batched decode only; O(N), fixed 2D; excludes calibration and physical resources'})
    save_csv(out/'table_1_all_results.csv',aggregate)
    save_csv(out/'runtime_samples.csv',[{'method':name,'repeat':i,'microseconds_per_trial':value} for name,values in times.items() for i,value in enumerate(values)])
    # Independent deterministic Gaussian integration of linear-decoder Fe.
    nodes,weights=np.polynomial.hermite.hermgauss(40)
    grid=np.array([(x,y) for x in nodes for y in nodes])*np.sqrt(2)
    quadrature_weights=np.outer(weights,weights).reshape(-1)/np.pi
    analytic=[]
    for j,sigma in enumerate(sigmas):
        s=noise_covariance(sigma,protocol['noise_axis_ratio'],protocol['noise_angle_radians'])
        for name in ('No recovery','ML','Geometry','Conventional posterior'):
            values=[]
            for k in range(seeds):
                if name=='ML':w=np.eye(2)
                elif name=='No recovery':w=np.zeros((2,2))
                else:w=np.linalg.solve(calibrated[k,j]+r,calibrated[k,j]).T
                c=(np.eye(2)-w)@s@(np.eye(2)-w).T+w@r@w.T
                eigen,vectors=np.linalg.eigh(c)
                gaussian=grid@(vectors*np.sqrt(np.maximum(eigen,0))).T
                values.append(float(quadrature_weights@metrics(gaussian)['entanglement_fidelity']))
            measured=next(x['entanglement_fidelity'] for x in aggregate if x['sigma']==sigma and x['method']==name)
            analytic.append({'sigma':sigma,'method':name,'quadrature_fe':float(np.mean(values)),
                             'monte_carlo_fe':measured,'mc_minus_quadrature':measured-float(np.mean(values))})
    save_csv(out/'analytic_reference.csv',analytic)
    physical=physical_validation(out)
    if max(audit_max.values())>1e-12:raise AssertionError('representation/equivalence check failed')
    npz_hash=hashlib.sha256((out/'raw_trials.npz').read_bytes()).hexdigest()
    record={"author":"Jason D. Dutton","date":"2026-10-06","protocol":protocol,
            "protocol_sha256":hashlib.sha256((ROOT/'protocol.json').read_bytes()).hexdigest(),
            "numpy_version":np.__version__,"python":sys.version,"platform":platform.platform(),
            "processor":platform.processor(),"logical_cpu_count":os.cpu_count(),
            "thread_environment":{v:os.environ.get(v) for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')},
            "held_out_trials":len(sigmas)*seeds*n,"calibration_injections":len(sigmas)*seeds*protocol['calibration_per_seed_per_noise'],
            "payload_mean_input_photons":0.5,"pilot_input_photons_per_trial":protocol['pilot_modes_per_trial']*abs(protocol['pilot_alpha'])**2,
            "controls":audit_max,"physical_validation":physical,"setup_seconds_all_decoders_by_noise":setup_times,
            "dataset_sha256":npz_hash,"wall_time_seconds":time.perf_counter()-started,
            "all_evaluation_uses_paired_trials":True,"decoder_access_to_test_truth":False,
            "geometry_is_conventional_gaussian_posterior":True,"no_hardware_or_fault_tolerance_claim":True}
    (out/'run_record.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(json.dumps({'trials':record['held_out_trials'],'controls':audit_max,'physical_validation':physical,'seconds':record['wall_time_seconds']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'results')
    args=parser.parse_args();run(json.loads((ROOT/'protocol.json').read_text()),args.output)
