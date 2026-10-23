"""Matched-information geometric recovery and oscillator convergence experiments."""
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import platform
import subprocess
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
spec=importlib.util.spec_from_file_location('repo_projection',ROOT/'SEFI-PY Evolution/qec/projection.py')
projection=importlib.util.module_from_spec(spec); spec.loader.exec_module(projection)
CONFIG=dict(seeds=[1701,1702,1703],samples_per_seed=10000,sigmas=[.001,.01,.1,.5],
            noise=['isotropic_gaussian','radial_gaussian','tangent_gaussian','student_t3'],
            frame_errors_radians=[0,.05,.2],duration=100.,steps=[.2,.1,.05,.025,.0125])

def project(points, radius=1., frame_error=0.):
    theta=np.arctan2(points[:,1],points[:,0])
    target=np.column_stack((radius*np.cos(theta),radius*np.sin(theta),np.zeros(len(theta))))
    normal=np.column_stack((np.cos(theta+frame_error),np.sin(theta+frame_error),np.zeros(len(theta))))
    b=np.tile([0.,0.,1.],(len(theta),1))
    delta=points-target
    dn=np.sum(delta*normal,axis=1); db=delta[:,2]
    return projection.project_to_manifold(points,normal,b,dn[:,None],db[:,None])

def conventional(points, radius=1.):
    theta=np.arctan2(points[:,1],points[:,0])
    return np.column_stack((radius*np.cos(theta),radius*np.sin(theta),np.zeros(len(theta))))

def oscillator(dt, method, duration=100.):
    n=round(duration/dt); dt=duration/n; x=np.array([1.,0.])
    def f(q): return np.array([q[1],-q[0]])
    t=time.perf_counter()
    for _ in range(n):
        if method in ['euler','projected_euler']:
            x=x+dt*f(x)
            if method=='projected_euler': x=x/np.linalg.norm(x)
        elif method=='midpoint':
            x=x+dt*f(x+dt/2*f(x))
        elif method=='rk4':
            a=f(x); b=f(x+dt*a/2); c=f(x+dt*b/2); d=f(x+dt*c)
            x=x+dt*(a+2*b+2*c+d)/6
    return dict(dt=dt,method=method,steps=n,seconds=time.perf_counter()-t,
                state_error=float(np.linalg.norm(x-[np.cos(duration),-np.sin(duration)])),
                energy_error=float(abs(x@x/2-.5)),
                phase_error=float(abs(np.angle(complex(x[0],-x[1])*np.exp(-1j*duration)))))

def main(out=None):
    out=out or Path(__file__).parent/'results'; out.mkdir(parents=True,exist_ok=True)
    rows=[]; raw={}
    for sigma in CONFIG['sigmas']:
        for noise in CONFIG['noise']:
            loss={}; timings={}
            for seed in CONFIG['seeds']:
                rng=np.random.default_rng(seed); n=CONFIG['samples_per_seed']
                theta=rng.uniform(-np.pi,np.pi,n)
                truth=np.column_stack((np.cos(theta),np.sin(theta),np.zeros(n)))
                tangent=np.column_stack((-np.sin(theta),np.cos(theta),np.zeros(n)))
                if noise=='radial_gaussian': perturb=sigma*rng.normal(size=(n,1))*truth
                elif noise=='tangent_gaussian': perturb=sigma*rng.normal(size=(n,1))*tangent
                elif noise=='student_t3': perturb=sigma*rng.standard_t(3,size=(n,3))/np.sqrt(3)
                else: perturb=sigma*rng.normal(size=(n,3))
                y=truth+perturb
                methods={'uncorrected':lambda:y,'nearest_circle':lambda:conventional(y),
                         'wrong_radius_1.05':lambda:project(y,1.05)}
                for angle in CONFIG['frame_errors_radians']:
                    methods[f'sefi_frame_{angle}']=lambda angle=angle:project(y,frame_error=angle)
                for name,fn in methods.items():
                    start=time.perf_counter(); estimate=fn(); seconds=time.perf_counter()-start
                    loss.setdefault(name,[]).append(np.sum((estimate-truth)**2,axis=1))
                    timings[name]=timings.get(name,0)+seconds
            base=np.concatenate(loss['nearest_circle'])
            for name,values in loss.items():
                values=np.stack(values); flat=values.ravel(); delta=flat-base
                rows.append(dict(sigma=sigma,noise=noise,method=name,samples=len(flat),
                    mse=float(flat.mean()),mse_mean_se=float(flat.std(ddof=1)/np.sqrt(len(flat))),
                    paired_mse_delta_vs_nearest=float(delta.mean()),
                    paired_delta_se=float(delta.std(ddof=1)/np.sqrt(len(flat))),
                    seed_mses=values.mean(axis=1).tolist(),seconds=timings[name]))
                raw[f'{sigma}_{noise}_{name}']=values
    conv=[oscillator(dt,method) for dt in CONFIG['steps'] for method in ['euler','projected_euler','midpoint','rk4']]
    # Existing formula assumes an orthonormal frame. Compare a general plane
    # projection that uses the Gram matrix for a nonorthogonal frame.
    p=np.array([.4,.8,1.2]); n=np.array([1.,0.,0.]); b=np.array([.6,.8,0.])
    A=np.column_stack((n,b)); d=A.T@p
    naive=projection.project_to_manifold(p,n,b,*d)
    gram=p-A@np.linalg.solve(A.T@A,d)
    diagnostic=dict(nonorthogonal_frame_repo_residual=float(np.linalg.norm(A.T@naive)),
                    gram_projection_residual=float(np.linalg.norm(A.T@gram)),
                    frame_condition_number=float(np.linalg.cond(A)))
    for name,obj in [('config',CONFIG),('results',rows),('convergence',conv),('diagnostics',diagnostic)]:
        (out/f'{name}.json').write_text(json.dumps(obj,indent=2))
    np.savez_compressed(out/'paired_squared_errors.npz',**raw)
    meta=dict(python=sys.version,numpy=np.__version__,platform=platform.platform(),
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        projection_sha256=hashlib.sha256((ROOT/'SEFI-PY Evolution/qec/projection.py').read_bytes()).hexdigest())
    (out/'environment.json').write_text(json.dumps(meta,indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for method in ['uncorrected','nearest_circle','sefi_frame_0','wrong_radius_1.05']:
        r=[x for x in rows if x['noise']=='isotropic_gaussian' and x['method']==method]
        axes[0].loglog([x['sigma'] for x in r],[x['mse'] for x in r],'o-',label=method)
    axes[0].set(xlabel='Noise standard deviation',ylabel='Recovery MSE',title='Known unit circle; matched observations')
    for method in ['euler','projected_euler','midpoint','rk4']:
        r=[x for x in conv if x['method']==method]
        axes[1].loglog([x['dt'] for x in r],[x['state_error'] for x in r],'o-',label=method)
    axes[1].set(xlabel='Time step',ylabel='State error at T=100',title='Dimensionless harmonic oscillator')
    for ax in axes: ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(out/'recovery_and_convergence.png',dpi=160); plt.close(fig)
    print(json.dumps(dict(cells=len(rows),convergence_runs=len(conv),diagnostic=diagnostic),indent=2))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path); args=parser.parse_args()
    main(args.output)
