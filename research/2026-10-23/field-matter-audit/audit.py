"""Action audit and charge continuation. Existing solvers/data are read-only.

This is a classical two-complex-scalar model; no particle identification.
"""
import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import time
import numpy as np
import scipy
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SOURCE = REPO/'research/2026-10-19/vortex-extended'
spec = importlib.util.spec_from_file_location('preserved_radial_model', SOURCE/'vacuum_spherical.py')
legacy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(legacy)
Radial = legacy.Radial


def potential(rho, eta, lam):
    return rho*rho/4 + rho/2 + 2*rho*eta + lam*(eta-3/lam)**2/4


def minimized_ratio(rho, lam):
    eta = np.maximum((3-4*np.asarray(rho))/lam, 0)
    return 2*potential(rho, eta, lam)/rho


def exact_square_identity(lam):
    """Exact polynomial coefficients, not a floating-point identity test.

    U - m_phi^2*rho/2 = lam/4*(eta-3/lam+4*rho/lam)^2
                              + (1/4-4/lam)*rho^2.
    Dictionaries index monomials by (power of rho, power of eta).
    """
    lam = Fraction(lam)
    left = {(2, 0): Fraction(1, 4), (1, 1): Fraction(2),
            (0, 2): lam/4, (0, 1): Fraction(-3, 2),
            (0, 0): Fraction(9, 4)/lam, (1, 0): -6/lam}
    a, b, c = Fraction(4)/lam, Fraction(1), -3/lam
    right = {(2, 0): lam*a*a/4 + Fraction(1, 4)-4/lam,
             (1, 1): lam*2*a*b/4, (0, 2): lam*b*b/4,
             (0, 1): lam*2*b*c/4, (0, 0): lam*c*c/4,
             (1, 0): lam*2*a*c/4}
    return left == right


def analytic_audit():
    records = []
    for lam in (4., 12., 15., 16., 25.):
        numerical = minimize_scalar(lambda log_rho: minimized_ratio(np.exp(log_rho), lam),
                                    bounds=(-24, 6), method='bounded')
        mass2 = 1+12/lam
        infimum = 1+3/np.sqrt(lam) if lam < 16 else mass2
        records.append({'lambda': lam, 'm_phi_squared': mass2,
                        'infimum_ratio': float(infimum), 'numeric_infimum': float(numerical.fun),
                        'subthreshold_well': lam < 16,
                        'exact_square_identity': exact_square_identity(int(lam))})
    return {'scope': 'Independent recheck of already documented potential bound; not a new theorem',
            'vacuum_real_component_mass_squared': {'phi_real': 2., 'phi_imag': 2.,
                'sigma_radial': 6., 'sigma_phase': 0.},
            'mass_spectrum_lambda': 12., 'coupling_audit': records}


def solve_profile(Q, h, L=80., tolerance=1e-9, maximum=12000):
    model = Radial(Q=float(Q), h=float(h), L=float(L), lam=12.)
    # A supplied profile is an initial guess, not a presumed solution at Q.
    data = np.load(SOURCE/'vacuum_spherical_h0p1_L80.npz', allow_pickle=False)
    scaled_r = model.r/(Q/100000.)**(1/3)
    f = np.interp(scaled_r, data['r'], data['f'], left=data['f'][0], right=0.)
    s = np.interp(scaled_r, data['r'], data['s'], left=data['s'][0], right=model.vac)
    energy = model.energy(f, s)
    start = time.perf_counter()
    reason = 'iteration_limit'
    for iteration in range(maximum):
        gf, gs, omega = model.gradient(f, s)
        residual = float(max(abs(gf).max(), abs(gs).max()))
        if residual < tolerance:
            reason = 'residual_tolerance'; break
        df = model.precondition(gf, 1+3*f*f+4*s*s+omega*omega)
        ds = model.precondition(gs, 1+abs(-3+4*f*f+3*model.lam*s*s))
        descent = float(np.sum(model.V*(gf*df+gs*ds)))
        step = 1.
        while step > 2**-24:
            fn, sn = f+step*df, s+step*ds
            candidate = model.energy(fn, sn)
            # At machine resolution, allow equality but never energy increase.
            if candidate <= energy-1e-4*step*descent:
                break
            step /= 2
        if step <= 2**-24:
            reason = 'line_search_stalled'; break
        f, s, energy = fn, sn, candidate
    gf, gs, omega = model.gradient(f, s)
    residual = float(max(abs(gf).max(), abs(gs).max()))
    I = np.sum(model.V*f*f)
    C = Q*Q/(2*I)
    U = float(np.sum(model.V*potential(f*f, s*s, model.lam)))
    K = energy-U-C
    result = {'Q': float(Q), 'h': float(h), 'L': float(L), 'lambda': model.lam,
              'energy': energy, 'energy_per_charge': energy/Q, 'omega': omega,
              'residual': residual, 'reason': reason, 'iterations': iteration,
              'requested_tolerance': tolerance, 'converged': residual < tolerance,
              'usable_stationarity_1e-7': residual < 1e-7,
              'K': float(K), 'U': U, 'C': float(C),
              'relative_virial_residual': float((K+3*U-3*C)/energy),
              'seconds': time.perf_counter()-start}
    return result, model, f, s


def branch_summary(records):
    summaries = []
    for h in sorted(set(r['h'] for r in records)):
        rows = sorted((r for r in records if r['h'] == h), key=lambda r: r['Q'])
        by_Q = {r['Q']: r for r in rows}
        slopes = []
        for Q in (50000., 100000., 150000.):
            if Q-1000 not in by_Q or Q+1000 not in by_Q:
                continue
            center, lower, upper = by_Q[Q], by_Q[Q-1000], by_Q[Q+1000]
            dEdQ = (upper['energy']-lower['energy'])/2000
            dw = (upper['omega']-lower['omega'])/2000
            slopes.append({'Q': Q, 'central_dE_dQ': dEdQ, 'omega': center['omega'],
                           'relative_identity_error': abs(dEdQ-center['omega'])/abs(center['omega']),
                           'domega_dQ': dw,
                           'eligible': all(x['usable_stationarity_1e-7'] for x in (lower, center, upper))})
        fission = []
        for parent, a, b in ((50000., 25000., 25000.), (100000., 50000., 50000.),
                             (150000., 50000., 100000.), (200000., 100000., 100000.)):
            if all(q in by_Q for q in (parent, a, b)):
                gap = by_Q[a]['energy']+by_Q[b]['energy']-by_Q[parent]['energy']
                fission.append({'parent_Q': parent, 'daughter_Qs': [a, b],
                                'daughter_energy_minus_parent': gap,
                                'fraction_of_parent_energy': gap/by_Q[parent]['energy'],
                                'eligible': all(by_Q[q]['usable_stationarity_1e-7'] for q in (parent, a, b)),
                                'scope': 'Selected same-branch separated candidates, not all fission channels'})
        summaries.append({'h': h, 'charge_derivatives': slopes, 'selected_fission': fission})
    return summaries


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    provenance = {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
                  'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
                      (SOURCE/'vacuum_spherical.py', SOURCE/'vacuum_spherical_h0p1_L80.npz', Path(__file__))},
                  'scope': 'Classical radial lambda=12 continuation. No new physics or observed-particle identification.'}
    (args.output/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    (args.output/'analytic.json').write_text(json.dumps(analytic_audit(), indent=2)+'\n')
    records = []
    charges = [25000., 49000., 50000., 51000., 75000., 99000., 100000., 101000.,
               125000., 149000., 150000., 151000., 200000.]
    for h in (.4, .2):
        for Q in charges:
            result, model, f, s = solve_profile(Q, h)
            records.append(result)
            np.savez_compressed(args.output/f'Q{int(Q)}_h{h}.npz', r=model.r, f=f, s=s)
            (args.output/'profiles.json').write_text(json.dumps(records, indent=2)+'\n')
            print(json.dumps(result), flush=True)
    (args.output/'branch_summary.json').write_text(json.dumps(branch_summary(records), indent=2)+'\n')


if __name__ == '__main__':
    main()
