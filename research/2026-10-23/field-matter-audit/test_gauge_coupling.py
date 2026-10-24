from fractions import Fraction
import numpy as np
from gauge_coupling import findings, gauge_test, covariant_velocity, covariant_charge


def test_interior_source_gauge_obstruction():
    result = findings()
    assert abs(result['constant_charge']['direct_action_change']) < 1e-12
    assert abs(result['changing_charge']['direct_action_change']-float(Fraction(28, 75))) < 1e-12
    assert result['changing_charge']['endpoint_contribution'] == 0
    assert result['changing_charge']['identity_error'] < 1e-12


def test_open_segment_endpoint_variation():
    curve = lambda t: np.array([t, 0., 0., 0.])
    tangent = lambda t: np.array([1., 0., 0., 0.])
    chi = lambda x: x[0]**2
    gradient = lambda x: np.array([2*x[0], 0., 0., 0.])
    result = gauge_test(curve, tangent, lambda t: 2., lambda t: 0., lambda x: np.zeros(4),
                       chi, gradient, 0., 2., .7)
    assert abs(result['direct_action_change']+5.6) < 1e-12
    assert result['interior_contribution'] == 0
    assert result['identity_error'] < 1e-12


def test_local_covariance_invariant_charge_and_required_quadratic_term():
    rng = np.random.default_rng(164)
    for _ in range(30):
        z = complex(*rng.normal(size=2)); velocity = complex(*rng.normal(size=2))
        A, chi, chi_dot = rng.normal(size=3); g = .7
        phase = np.exp(1j*g*chi)
        transformed_z = phase*z
        transformed_velocity = phase*(velocity+1j*g*chi_dot*z)
        D = covariant_velocity(z, velocity, A, g)
        transformed_D = covariant_velocity(transformed_z, transformed_velocity, A+chi_dot, g)
        assert abs(transformed_D-phase*D) < 1e-12
        assert abs(covariant_charge(z, velocity, A, g)-covariant_charge(transformed_z, transformed_velocity, A+chi_dot, g)) < 1e-12
        ordinary_charge = np.imag(np.conj(z)*velocity)
        expanded = .5*abs(velocity)**2-g*A*ordinary_charge+.5*g*g*A*A*abs(z)**2
        assert abs(.5*abs(D)**2-expanded) < 1e-12
        eps = 1e-5
        derivative = (.5*abs(covariant_velocity(z, velocity, A+eps, g))**2-
                      .5*abs(covariant_velocity(z, velocity, A-eps, g))**2)/(2*eps)
        assert abs(derivative+g*covariant_charge(z, velocity, A, g)) < 1e-9
