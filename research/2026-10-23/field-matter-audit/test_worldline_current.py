from fractions import Fraction
import numpy as np
from worldline_current import weak_divergence, source_action, weak_current, endpoint_balance, findings


def test_endpoint_signs_against_analytic_polynomial():
    curve = lambda t: np.array([t, .2*t, .1*t, 0.])
    tangent = lambda t: np.array([1., .2, .1, 0.])
    test = lambda x: x[0]**2+x[1]*x[2]
    gradient = lambda x: np.array([2*x[0], x[2], x[1], 0.])
    actual, _ = weak_divergence(curve, tangent, lambda t: 2., gradient, 0., 2.)
    expected = source_action(curve, lambda t: 2., lambda t: 0., test, 0., 2.)
    assert abs(actual+8.16) < 1e-12
    assert expected['initial_endpoint'] == 0
    assert abs(expected['final_endpoint']+8.16) < 1e-12
    assert abs(actual-expected['total']) < 1e-12


def test_interior_torque_is_not_conserved_charge():
    records = findings()['examples']
    assert abs(records[0]['weak_divergence']) < 1e-12
    assert abs(records[1]['weak_divergence']-float(Fraction(8, 15))) < 1e-12
    assert all(row['identity_absolute_error'] < 1e-12 for row in records)


def test_reparametrization_current_and_divergence():
    curve = lambda t: np.array([t, .2*np.sin(t), .1*np.cos(t), 0.])
    tangent = lambda t: np.array([1., .2*np.cos(t), -.1*np.sin(t), 0.])
    charge = lambda t: 1+.1*t*t
    covector = lambda x: np.array([x[0]**2+x[1], np.exp(.1*x[0]), x[0]*x[2], 0.])
    gradient = lambda x: np.array([2*x[0], x[2], x[1], 0.])
    g = lambda u: u+.2*np.sin(u)
    gp = lambda u: 1+.2*np.cos(u)
    for test_function in (covector, gradient):
        reference, _ = weak_current(curve, tangent, charge, test_function, g(-1.), g(1.))
        transformed, _ = weak_current(lambda u: curve(g(u)), lambda u: tangent(g(u))*gp(u),
                                      lambda u: charge(g(u)), test_function, -1., 1.)
        assert abs(reference-transformed) < 1e-11


def test_vertex_balance_and_pair_endpoint_cancellation():
    start, vertex = (0, 0, 0, 0), (1, 0, 0, 0)
    a, b = (2, Fraction(1, 5), 0, 0), (2, Fraction(-1, 5), 0, 0)
    balanced = endpoint_balance([(start, vertex, 3), (vertex, a, 1), (vertex, b, 2)])
    unbalanced = endpoint_balance([(start, vertex, 3), (vertex, a, 1), (vertex, b, 1)])
    assert balanced[vertex] == 0
    assert unbalanced[vertex] == -1
    paired = endpoint_balance([(vertex, a, 1), (vertex, b, -1)])
    assert paired[vertex] == 0
    assert paired[a] == -1 and paired[b] == 1
