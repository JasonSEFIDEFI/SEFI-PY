"""Full radial amplitude/phase linearization in Phi's rotating frame.

Fixed omega here labels the representation, not fixed-charge relaxation.
No finite spectrum is promoted to a continuum nonlinear stability theorem.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.linalg import eig, eigh
from audit import Radial, solve_profile


def matrices(model, f, s):
    n = len(f); root = np.sqrt(model.V)
    omega = model.Q/np.sum(model.V*f*f)
    kin = np.diag(model.diag/model.V)
    off = -model.links/(root[:-1]*root[1:])
    kin += np.diag(off, 1)+np.diag(off, -1)
    plus = kin+np.diag(1+3*f*f+4*s*s-omega**2)
    minus = kin+np.diag(1+f*f+4*s*s-omega**2)
    sigma = kin+np.diag(-3+4*f*f+3*model.lam*s*s)
    sigma_phase = kin+np.diag(-3+4*f*f+model.lam*s*s)
    cross = np.diag(8*f*s); zero = np.zeros((n, n)); eye = np.eye(n)
    K = np.block([[plus, zero, cross], [zero, minus, zero], [cross, zero, sigma]])
    G = np.block([[zero, -2*omega*eye, zero], [2*omega*eye, zero, zero], [zero, zero, zero]])
    A = np.block([[np.zeros_like(K), np.eye(3*n)], [-K, -G]])
    local_real = np.block([[plus, cross], [cross, sigma]])
    u = np.r_[root*f, np.zeros(n)]
    constrained = local_real+4*omega**2/np.sum(model.V*f*f)*np.outer(u, u)
    return A, local_real, constrained, sigma_phase


def calculate(h):
    profile, model, f, s = solve_profile(100000., h)
    A, local, constrained, sigma_phase = matrices(model, f, s)
    start = time.perf_counter()
    values, vectors = eig(A, check_finite=False)
    residuals = np.linalg.norm(A@vectors-vectors*values, axis=0)
    scale = np.linalg.norm(A, ord=np.inf)+np.abs(values)
    normalized = residuals/scale
    # Pure Phi global phase mode, expressed in volume-weighted coordinates.
    phase = np.zeros(A.shape[0]); phase[len(f):2*len(f)] = np.sqrt(model.V)*f
    phase /= np.linalg.norm(phase)
    positive = np.where(values.real > 1e-7)[0]
    selected = sorted(positive, key=lambda i: -values[i].real)
    return {'h': h, 'Q': model.Q, 'L': model.L, 'dimension': A.shape[0],
            'profile': profile, 'seconds': time.perf_counter()-start,
            'real_grand_canonical_lowest': eigh(local, subset_by_index=(0, 2), eigvals_only=True).tolist(),
            'real_fixed_charge_lowest': eigh(constrained, subset_by_index=(0, 2), eigvals_only=True).tolist(),
            'sigma_phase_lowest': eigh(sigma_phase, subset_by_index=(0, 2), eigvals_only=True).tolist(),
            'max_real_part': float(values.real.max()), 'max_normalized_eigenpair_residual': float(normalized.max()),
            'positive_real_roots_above_1e-7': [{'real': float(values[i].real), 'imag': float(values[i].imag),
                'phase_overlap': float(abs(phase @ vectors[:, i])),
                'absolute_residual': float(residuals[i]), 'normalized_residual': float(normalized[i])} for i in selected],
            'near_zero_roots': [[float(z.real), float(z.imag)] for z in values[np.argsort(abs(values))[:6]]],
            'scope': 'All finite-grid radial coupled roots; decoupled Sigma phase tested by its curvature. Positive roots require symmetry/refinement audit.'}


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    results = []
    for h in (.8, .4):
        result = calculate(h); results.append(result)
        print(json.dumps(result), flush=True)
    args.output.write_text(json.dumps(results, indent=2)+'\n')


if __name__ == '__main__':
    main()
