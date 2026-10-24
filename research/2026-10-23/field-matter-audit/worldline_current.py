"""Weak distributional worldline-current checks, not bulk field dynamics."""
from fractions import Fraction
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad


def integrate(function, a, b):
    if not a < b:
        raise ValueError('ordered finite integration interval required')
    value, error = quad(function, a, b, epsabs=1e-11, epsrel=1e-11, limit=200)
    return float(value), float(error)


def weak_divergence(curve, tangent, charge, gradient, a, b):
    # <partial_mu j^mu, f> = - integral J dy^mu/du partial_mu f du.
    return integrate(lambda u: -charge(u)*np.dot(tangent(u), gradient(curve(u))), a, b)


def source_action(curve, charge, charge_derivative, test, a, b):
    bulk, error = integrate(lambda u: charge_derivative(u)*test(curve(u)), a, b)
    initial = charge(a)*test(curve(a))
    final = -charge(b)*test(curve(b))
    return {'initial_endpoint': float(initial), 'final_endpoint': float(final),
            'interior_source': bulk, 'total': float(initial+final+bulk), 'quadrature_error_estimate': error}


def weak_current(curve, tangent, charge, covector_test, a, b):
    return integrate(lambda u: charge(u)*np.dot(tangent(u), covector_test(curve(u))), a, b)


def endpoint_balance(segments):
    """Exact endpoint coefficients: outgoing charge minus incoming charge."""
    values = {}
    for start, end, charge in segments:
        for event, coefficient in ((tuple(start), Fraction(charge)), (tuple(end), -Fraction(charge))):
            values[event] = values.get(event, Fraction(0))+coefficient
    return values


def findings():
    curve = lambda t: np.array([t, .2*t, .1*t, 0.])
    tangent = lambda t: np.array([1., .2, .1, 0.])
    test = lambda x: x[0]**2*(2-x[0])**2
    gradient = lambda x: np.array([4*x[0]*(2-x[0])*(1-x[0]), 0., 0., 0.])
    examples = []
    for name, charge, derivative in [('constant_charge', lambda t: 2., lambda t: 0.),
                                      ('changing_charge', lambda t: 1+.5*t, lambda t: .5)]:
        actual, error = weak_divergence(curve, tangent, charge, gradient, 0., 2.)
        expected = source_action(curve, charge, derivative, test, 0., 2.)
        examples.append({'name': name, 'weak_divergence': actual, 'quadrature_error_estimate': error,
                         'source_decomposition': expected, 'identity_absolute_error': abs(actual-expected['total'])})
    return {'scope': 'Finite future-timelike inertial example; tests vanish at endpoints; not a field-generation simulation',
            'curve': 'y(t)=(t,0.2t,0.1t,0), 0<=t<=2',
            'test_on_curve': 't^2(2-t)^2, smooth compact cutoff understood outside the finite curve',
            'examples': examples,
            'exact_expected_changing_charge_result': '8/15',
            'interpretation': 'Conservation requires zero interior torque and accounted endpoints; current support remains a curve.'}


if __name__ == '__main__':
    path = Path(__file__).with_name('worldline-current-results.json')
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(findings(), indent=2)+'\n')
