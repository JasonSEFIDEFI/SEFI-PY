import numpy as np
from internal_symmetry import symmetry_conditions, potential, acceleration, charges, torque


def test_normalization_reveals_weaker_required_equalities():
    result = symmetry_conditions([1, 4, 9, 16], [2, 8, 27, 48])
    assert result['ordinary_origin_centered_pair_rotations'] == [False, False]
    assert result['kinetically_normalized_pair_rotations'] == [True, True]
    assert symmetry_conditions([1]*4, [1, 2, 3, 4])['kinetically_normalized_pair_rotations'] == [False, False]


def test_normalized_potential_rotation():
    h = np.array([1., 4., 9., 16.]); k = np.array([2., 8., 27., 48.])
    q = np.array([1., 2., 3., 4.]); rotated = np.array([-2., 1., -4., 3.])
    assert potential(q, h, k, 2., 3., [0]*4) == potential(rotated, h, k, 2., 3., [0]*4)
    assert potential(q, h, k, 2., 3., [1, 0, 2, 0]) != potential(rotated, h, k, 2., 3., [1, 0, 2, 0])


def test_rhs_gradient_and_charge_torque_independent_differences():
    rng = np.random.default_rng(63)
    for _ in range(20):
        h = rng.uniform(.5, 3, 4); k = rng.uniform(.5, 4, 4)
        q = rng.normal(size=4); velocity = rng.normal(size=4); coefficients = rng.normal(size=4)
        acc = acceleration(q, h, k, 1.2, .8, coefficients)
        eps = 1e-5; gradient = []
        for e in np.eye(4):
            gradient.append((potential(q+eps*e, h, k, 1.2, .8, coefficients)-potential(q-eps*e, h, k, 1.2, .8, coefficients))/(2*eps))
        np.testing.assert_allclose(acc, -np.array(gradient), atol=1e-9, rtol=1e-8)
        derivative = (charges(q+eps*velocity, velocity+eps*acc)-charges(q-eps*velocity, velocity-eps*acc))/(2*eps)
        np.testing.assert_allclose(derivative, torque(q, h, k, 1.2, .8, coefficients), atol=1e-9, rtol=1e-8)


def test_shifted_symmetry_survives_constant_prescribed_geometry():
    h = np.array([1., 4., 9., 16.]); k = np.array([2., 8., 27., 48.])
    q = np.array([1., 2., 3., 4.]); coefficients = [1., 0., 2., 0.]
    omega2 = k/h
    source = np.array([4/np.sqrt(h[2]), 18/np.sqrt(h[3])])
    center = -source/omega2[2:]
    acc = acceleration(q, h, k, 2., 3., coefficients)
    shifted = q[2:]-center
    assert abs(shifted[0]*acc[3]-shifted[1]*acc[2]) < 1e-12
    assert abs(torque(q, h, k, 2., 3., coefficients)[1]) > 1
