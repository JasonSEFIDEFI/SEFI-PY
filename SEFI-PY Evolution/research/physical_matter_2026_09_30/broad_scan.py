"""Bound-informed N=0 falsification scan; no automatic instability acceptance.

Retains the original full-Z operator and checkpoints. Restartable per target.
"""
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path
import os

# Avoid BLAS oversubscription on sparse sequential jobs. Explicit user settings
# take precedence, and these variables must precede numerical-library imports.
for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(variable, '1')

import numpy as np
import scipy
from scipy import sparse as sp
from scipy.sparse.linalg import LinearOperator, eigs, eigsh, splu, ArpackNoConvergence
from full_z_spectral import load_profile, stationary, spectral_operator, measure_null


def weights(n, h, m):
    r = np.arange(0 if m == 0 else 1, n, dtype=float)
    if m == 0:
        r[0] = 1 / 8  # discrete axis detailed balance: w0*4 = w1/2
    return np.tile(np.repeat(r*h*h, 2*n-1), 4)


def structure_bound(b, g, w):
    root = np.sqrt(w)
    bh = sp.diags(root) @ b @ sp.diags(1/root)
    gh = sp.diags(root) @ g @ sp.diags(1/root)
    herm = float(abs(bh-bh.getH()).max())
    skew = float(abs(gh+gh.getH()).max())
    if max(herm, skew) > 1e-10:
        raise RuntimeError('Energy bound requires weighted Hermitian B and skew-Hermitian G')
    vals, vec = eigsh(bh, k=1, which='LA', tol=1e-10,
        v0=np.random.default_rng(20260930).standard_normal(bh.shape[0]))
    residual = float(np.linalg.norm(bh@vec[:, 0]-vals[0]*vec[:, 0]))
    return {'B_hermitian_defect': herm, 'G_skew_defect': skew,
            'lambda_max_B_estimate': float(vals[0]), 'bound_eigenpair_residual': residual,
            'bound_start_seed':20260930,
            'growth_disk_radius_estimate': float(np.sqrt(max(0, vals[0]))),
            'certified': False,
            'meaning': 'All nonimaginary eigenvalues obey |sigma|^2 <= lambda_max(B); numerical extremal estimate is not a certified upper bound.'}


def schur_inverse(b, g, shift):
    """Exact block elimination of (A-shift I), independent factorization path."""
    size = b.shape[0]
    q = b + shift*g - shift**2*sp.eye(size, format='csc')
    lu = splu(q.tocsc())
    def solve(v):
        f, z = v[:size], v[size:]
        y = lu.solve(z - g@f + shift*f)
        return np.r_[y, f+shift*y]
    return LinearOperator((2*size, 2*size), matvec=solve, dtype=np.complex128)


def diagnose(value, v, full, b, g, modes, w, n, h, m):
    size = b.shape[0]
    y, vel = v[:size], v[size:]
    norm = lambda x: float(np.sqrt(np.vdot(x, w*x).real))
    q = b@y + value*(g@y) - value**2*y
    denom = norm(b@y) + abs(value)*norm(g@y) + abs(value)**2*norm(y)
    energy = (w*abs(y)**2).reshape(4, n-(m != 0), 2*n-1).sum(axis=0)
    r = np.arange(0 if m == 0 else 1, n)*h
    z = np.arange(-n+1, n)*h
    boundary = (r[:, None] >= .9*n*h) | (abs(z)[None, :] >= .9*n*h)
    peak = np.unravel_index(np.argmax(energy), energy.shape)
    def overlap(a, v):
        return float(abs(np.vdot(a, w*v))/(norm(a)*norm(v)))
    overlaps = {}
    for name, tangent in modes.items():
        vv, tt = y.copy(), tangent.copy()
        if name == 'longitudinal_translation':
            phase = modes['carrier_phase']
            vv -= phase*np.vdot(phase, w*vv)/np.vdot(phase, w*phase)
            tt -= phase*np.vdot(phase, w*tt)/np.vdot(phase, w*phase)
        overlaps[name] = overlap(tt, vv)
    # Variation across neighboring nodes measures resolution; not a pass by itself.
    yy = y.reshape(4, n-(m != 0), 2*n-1)
    variation = sum(np.sum(abs(np.diff(yy, axis=ax))**2) for ax in (1, 2))/np.sum(abs(yy)**2)
    return {'real': float(value.real), 'imag': float(value.imag),
            'first_order_residual': float(np.linalg.norm(full@v-value*v)/np.linalg.norm(v)),
            'quadratic_absolute_residual': norm(q)/norm(y),
            'quadratic_scaled_residual': norm(q)/max(denom, 1e-30),
            'velocity_consistency': norm(vel-value*y)/norm(y),
            'boundary_norm_fraction': float(np.sqrt(energy[boundary].sum()/energy.sum())),
            'peak_r_z': [float(r[peak[0]]), float(z[peak[1]])],
            'rms_r_z': [float(np.sqrt((energy*r[:, None]**2).sum()/energy.sum())),
                          float(np.sqrt((energy*z[None, :]**2).sum()/energy.sum()))],
            'neighbor_variation': float(variation), 'weighted_symmetry_overlaps': overlaps,
            'candidate': bool(value.real > 1e-6), 'accepted_instability': False}


def target_grid(radius):
    # Off-axis targets catch oscillatory growth as well as real growth.
    # No finite k-nearest target grid is asserted complete.
    extent = max(.25, radius*1.1)
    return [1e-4+0j] + [x*extent+1j*y*extent
        for x in (.15, .6) for y in (-1., -.5, 0., .5, 1.)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('profile', type=Path)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--h', type=float, required=True)
    ap.add_argument('--m', type=int, nargs='+', default=[0, 1, 2, 3, 4])
    ap.add_argument('--k', type=int, default=12)
    ap.add_argument('--shifts', type=complex, nargs='+')
    ap.add_argument('--strategies', nargs='+', choices=['schur', 'full', 'LR'], default=['schur'])
    ap.add_argument('--maxiter', type=int, default=1200)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    source_hash = hashlib.sha256(args.profile.read_bytes()).hexdigest()
    previous_meta = args.output/'metadata.json'
    if previous_meta.exists():
        previous = json.loads(previous_meta.read_text())
        if (previous['source_sha256'],previous['n'],previous['h']) != (source_hash,args.n,args.h):
            raise ValueError('Existing output has a different source/grid; use a new directory')
    n,h,fields,c,nu,lam,N = load_profile(args.profile, args.n, args.h)
    if N != 0 or n < 4 or h <= 0 or any(m < 0 for m in args.m):
        raise ValueError('Only N=0, n>=4, h>0, m>=0 supported')
    fields, hist = stationary(fields,n,h,c,nu,lam,30)
    if hist[-1] > 1e-9 or not np.isfinite(hist[-1]) or fields[2].max() < 1e-4:
        raise RuntimeError('Stationary profile convergence/nontriviality gate failed')
    np.savez_compressed(args.output/'profile.npz', r=np.arange(n+1)*h,
        z=np.arange(-n,n+1)*h, u=fields[0], w=fields[1], s=fields[2], c=c, nu=nu, lambda_sigma=lam)
    meta = {'source_sha256': source_hash,
        'n':n, 'h':h, 'L':n*h, 'N':N, 'c':c, 'nu':nu, 'lambda_sigma':lam,
        'stationary_history':hist, 'python':platform.python_version(),
        'numpy':np.__version__, 'scipy':scipy.__version__,
        'acceptance':'No candidate accepted automatically; all six user gates required',
        'sectors':{}}
    for m in args.m:
        full,b,g,modes,sl = spectral_operator(fields,n,h,c,nu,lam,m)
        w = weights(n,h,m)
        bound = structure_bound(b,g,w)
        meta['sectors'][str(m)] = {'dimension': full.shape[0], 'bound':bound,
            'neutral_residuals':{name:measure_null(b,y) for name,y in modes.items()}}
        (args.output/'metadata.json').write_text(json.dumps(meta,indent=2))
        shifts = args.shifts or target_grid(bound['growth_disk_radius_estimate'])
        for strategy in args.strategies:
            for idx,shift in enumerate(shifts if strategy != 'LR' else [None]):
                stem = f'm{m}_{strategy}_{idx:02d}'
                dest = args.output/(stem+'.json')
                if dest.exists():
                    previous = json.loads(dest.read_text())
                    requested_shift = None if shift is None else [shift.real,shift.imag]
                    if previous['shift'] != requested_shift or previous['k_requested'] != args.k:
                        raise ValueError(f'{dest}: existing target/k differ; use a new output directory')
                    continue
                start = time.perf_counter()
                status, error = 'complete', None
                ev,vec = np.empty(0,complex), np.empty((full.shape[0],0),complex)
                try:
                    opts = dict(k=args.k,tol=1e-10,maxiter=args.maxiter,
                        ncv=max(2*args.k+1,40),
                        v0=np.random.default_rng(20260930+idx).standard_normal(full.shape[0]))
                    if strategy == 'LR':
                        ev,vec = eigs(full,which='LR',**opts)
                    else:
                        inv = {'OPinv':schur_inverse(b,g,shift)} if strategy == 'schur' else {}
                        ev,vec = eigs(full,sigma=shift,which='LM',**inv,**opts)
                except ArpackNoConvergence as exc:
                    status,error = 'partial', str(exc)
                    ev,vec = exc.eigenvalues,exc.eigenvectors
                except Exception as exc:
                    status,error = 'failed', repr(exc)
                records = [diagnose(value,v,full,b,g,modes,w,n,h,m) for value,v in zip(ev,vec.T)]
                # Retain representative vectors and every nonsymmetry growth
                # candidate; JSON retains diagnostics for every returned pair.
                keep = [j for j,rec in enumerate(records) if idx == 0 or
                    (rec['candidate'] and max(rec['weighted_symmetry_overlaps'].values(),default=0) < .9)]
                np.savez_compressed(args.output/(stem+'.npz'),eigenvalues=ev,
                    saved_indices=np.array(keep,dtype=int),eigenvectors=vec[:,keep])
                report = {'m':m,'strategy':strategy,'shift':None if shift is None else [shift.real,shift.imag],
                    'k_requested':args.k,'maxiter':args.maxiter,'ncv':max(2*args.k+1,40),
                    'tol':1e-10,'seed':20260930+idx,
                    'status':status,'error':error,'seconds':time.perf_counter()-start,
                    'eigenvalues':records}
                dest.write_text(json.dumps(report,indent=2))
                print(json.dumps({'run':str(args.output),'target':stem,'status':status,
                    'seconds':report['seconds'],'max_real':max([e.real for e in ev],default=None)}),flush=True)


if __name__ == '__main__':
    main()
