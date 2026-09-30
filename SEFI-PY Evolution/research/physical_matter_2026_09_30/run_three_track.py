"""Reproduce continuum diagnostics and independent centered continuation."""
import os
for variable in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(variable,'1')
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import scipy
from full_z_spectral import load_profile,stationary
from branch_audit import diagnostics


def solve(source,n,h,c,nu,lam,out,fields=None):
    out.mkdir(parents=True,exist_ok=True)
    if (out/'diagnostics.json').exists():
        row=json.loads((out/'diagnostics.json').read_text())
        if (row['n'],row['h'],row['c'],row['nu'])!=(n,h,c,nu):raise ValueError('Restart settings mismatch')
        with np.load(out/'profile.npz') as d:f=[d[k] for k in ('u','w','s')]
        return f,row
    if fields is None:_,_,fields,_,_,_,N=load_profile(source,n,h)
    else:fields=[f.copy() for f in fields];N=0
    start=time.perf_counter()
    row={'n':n,'h':h,'c':c,'nu':nu,'lambda':lam,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()}
    try:
        if N!=0:raise ValueError('Only N=0 supported')
        fields,hist=stationary(fields,n,h,c,nu,lam,35)
        row.update(newton_history=hist,stationary_residual=hist[-1])
        if not np.isfinite(hist[-1]) or hist[-1]>1e-9 or fields[2].max()<1e-4:raise RuntimeError('Stationarity/nontriviality gate failed')
        row.update(diagnostics(fields,n,h,c,nu,lam),status='converged')
        np.savez_compressed(out/'profile.npz',r=np.arange(n+1)*h,z=np.arange(-n,n+1)*h,u=fields[0],w=fields[1],s=fields[2],c=c,nu=nu,lambda_sigma=lam)
    except Exception as error:
        row.update(status='failed',error=repr(error))
        (out/'failure.json').write_text(json.dumps(row,indent=2));print(json.dumps(row),flush=True)
        return None,row
    row['seconds']=time.perf_counter()-start
    (out/'diagnostics.json').write_text(json.dumps(row,indent=2))
    print(json.dumps({k:row[k] for k in ('L','h','c','nu','status','seconds','E_excess','P_physical_z','Q_sigma','radius_carrier_rms','pohozaev_transverse_relative','pohozaev_longitudinal_relative')}),flush=True)
    return fields,row


def centered_checks(rows,center,parameter,steps):
    out=[];omega=np.sqrt(2)
    for delta in steps:
        low=rows.get(round(center[parameter]-delta,8));high=rows.get(round(center[parameter]+delta,8))
        if low is None or high is None or low['status']!='converged' or high['status']!='converged':continue
        derivative=lambda key:(high[key]-low[key])/(2*delta)
        terms=[derivative('E_excess'),center['v']*derivative('P_physical_z'),omega*derivative('Delta_Q_phi'),center['nu']/center['gamma']*derivative('Q_sigma')]
        first=terms[0]-sum(terms[1:])
        correction=-center['v']/(2*omega*center['gamma']**2)*center['pohozaev_longitudinal'] if parameter=='c' else 0.
        expected=.5*center['P_reduced_discrete_variational'] if parameter=='c' else -center['Q_sigma_discrete']
        expected_cont=.5*center['P_reduced_bare'] if parameter=='c' else -center['Q_sigma']
        out.append({'parameter':parameter,'step':delta,'first_law_terms':terms,'first_law_residual':first,
            'first_law_relative':first/max(sum(map(abs,terms)),1e-30),
            'finite_box_first_law_prediction':correction,'first_law_minus_box_prediction':first-correction,
            'discrete_envelope_derivative':derivative('G_discrete'),'discrete_envelope_expected':expected,
            'discrete_envelope_error':derivative('G_discrete')-expected,
            'continuum_quadrature_envelope_error':derivative('G_continuum_quadrature')-expected_cont,
            'dE_dP_raw':derivative('E_excess')/derivative('P_physical_z') if abs(derivative('P_physical_z'))>1e-12 else None})
    return out


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=Path('three_track_runs'))
    ap.add_argument('--track',choices=['convergence','continuation','all'],default='all')
    ap.add_argument('--fine-only',action='store_true');a=ap.parse_args()
    data=Path(__file__).resolve().parents[2]/'sim/vortex_console/data'
    source=data/'charged_ring_refined.npz';wide=data/'charged_ring_wide.npz';a.output.mkdir(parents=True,exist_ok=True)
    rows=[];checks=[]
    if a.track in ('convergence','all'):
        jobs=[(32,.5,source),(40,.4,source),(48,1/3,source),(48,.5,source),(60,.4,source),(72,1/3,source),(64,.5,wide),(80,.4,wide),(96,1/3,wide),(80,.5,wide),(100,.4,wide),(120,1/3,wide)]
        for n,h,p in jobs:
            _,row=solve(p,n,h,.6,.3,25,a.output/f'convergence_L{n*h:g}_h{h:.6f}');rows.append(row)
            (a.output/'convergence.json').write_text(json.dumps(rows,indent=2))
    if a.track in ('continuation','all'):
        grids=[(72,1/3)] if a.fine_only else [(32,.5),(48,.5),(60,.4)]
        for n,h in grids:
            prefix=a.output/f'continuation_L{n*h:g}_h{h:g}'
            seed,center=solve(source,n,h,.6,.3,25,prefix/'center')
            if seed is None:continue
            allrows=[center]
            choices=[('c',[-.1,-.08,-.06,-.04,-.02,-.01,.01,.02,.04,.06,.08,.1]),('nu',[-.04,-.02,-.01,.01,.02,.04])]
            if a.fine_only:choices=[(p,[-.02,-.01,-.005,-.0025,.0025,.005,.01,.02]) for p in ('c','nu')]
            for parameter,offsets in choices:
                branch={round(center[parameter],8):center}
                # Continue outward separately in each direction from the same center.
                for sign in (-1,1):
                    previous=seed
                    for offset in sorted([x for x in offsets if x*sign>0],key=abs):
                        c=round(.6+(offset if parameter=='c' else 0),8);nu=round(.3+(offset if parameter=='nu' else 0),8)
                        f,row=solve(source,n,h,c,nu,25,prefix/f'{parameter}_{c if parameter=="c" else nu:.4f}',previous)
                        branch[round(row[parameter],8)]=row;allrows.append(row)
                        if f is None:break
                        previous=f
                local=centered_checks(branch,center,parameter,[.02,.01,.005,.0025] if a.fine_only else [.02,.01])
                checks.extend([dict(x,L=n*h,h=h) for x in local])
                (a.output/('first_law_fine_checks.json' if a.fine_only else 'first_law_checks.json')).write_text(json.dumps(checks,indent=2))
                (prefix/'branch.json').write_text(json.dumps(allrows,indent=2))
            keys=['c','nu','L','h','E_excess','P_physical_z','Delta_Q_phi','Q_sigma','radius_carrier_rms','radius_carrier_peak','stationary_residual','status']
            with (prefix/'branch.csv').open('w',newline='') as handle:
                writer=csv.DictWriter(handle,fieldnames=keys);writer.writeheader();writer.writerows({k:r.get(k) for k in keys} for r in allrows)
    (a.output/'environment.json').write_text(json.dumps({'numpy':np.__version__,'scipy':scipy.__version__,'threads':{k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')}},indent=2))


if __name__=='__main__':main()
