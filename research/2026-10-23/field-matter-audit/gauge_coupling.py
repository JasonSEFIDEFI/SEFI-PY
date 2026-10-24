"""Consistency tests for explicitly proposed gauge couplings; not new physics."""
import json
from pathlib import Path
import numpy as np
from worldline_current import weak_current, source_action


def interaction(curve, tangent, charge, potential, a, b, g):
    value, error = weak_current(curve, tangent, charge, potential, a, b)
    return -g*value, abs(g)*error


def gauge_test(curve, tangent, charge, derivative, potential, chi, gradient, a, b, g):
    original, _ = interaction(curve, tangent, charge, potential, a, b, g)
    transformed, error = interaction(curve, tangent, charge, lambda x: potential(x)+gradient(x), a, b, g)
    source = source_action(curve, charge, derivative, chi, a, b)
    return {'direct_action_change': transformed-original, 'predicted_action_change': g*source['total'],
            'identity_error': abs(transformed-original-g*source['total']),
            'endpoint_contribution': g*(source['initial_endpoint']+source['final_endpoint']),
            'interior_contribution': g*source['interior_source'], 'quadrature_error_estimate': error}


def covariant_velocity(z, velocity, A_pullback, g):
    return velocity-1j*g*A_pullback*z


def covariant_charge(z, velocity, A_pullback, g):
    return float(np.imag(np.conj(z)*covariant_velocity(z, velocity, A_pullback, g)))


def findings():
    curve = lambda t: np.array([t, .2*t, .1*t, 0.])
    tangent = lambda t: np.array([1., .2, .1, 0.])
    potential = lambda x: np.array([.2*x[0], .1*x[1], -.3*x[2], 0.])
    chi = lambda x: x[0]**2*(2-x[0])**2
    gradient = lambda x: np.array([4*x[0]*(2-x[0])*(1-x[0]), 0., 0., 0.])
    return {'scope': 'Proposed S_int=-g integral j.A, g=0.7; chi vanishes at endpoints',
            'constant_charge': gauge_test(curve, tangent, lambda t: 2., lambda t: 0., potential, chi, gradient, 0., 2., .7),
            'changing_charge': gauge_test(curve, tangent, lambda t: 1+.5*t, lambda t: .5, potential, chi, gradient, 0., 2., .7),
            'exact_expected_changing_charge_variation': '28/75',
            'interpretation': 'Endpoint phases alone cannot repair a nonzero interior torque. Gauging an invariant internal sector changes its derivative and source current.'}


if __name__ == '__main__':
    path = Path(__file__).with_name('gauge-coupling-results.json')
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(findings(), indent=2)+'\n')
