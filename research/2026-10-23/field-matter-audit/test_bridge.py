from fractions import Fraction as F
from bridge_check import invariant_monomials, portal_identity, internal_potential


def test_complete_polynomial_basis():
    assert set(invariant_monomials()) == {(0, 0, 0, 0), (1, 1, 0, 0), (0, 0, 1, 1),
                                          (2, 2, 0, 0), (1, 1, 1, 1), (0, 0, 2, 2)}


def test_generic_identity_and_symmetry_cannot_choose_binding():
    for a, b, c, g in ((1, 12, 1, 4), (1, 25, 1, 4), (2, 9, 3, 2)):
        assert portal_identity(a, b, c, g, F(3, b))['identity_exact']
    assert F(portal_identity(1, 12, 1, 4, F(1, 4))['quartic_discriminant_ab_minus_g2']) < 0
    assert F(portal_identity(1, 25, 1, 4, F(3, 25))['quartic_discriminant_ab_minus_g2']) > 0


def test_naive_pairing_fails_for_allowed_coupling():
    assert internal_potential((0, 0, 1, 0)) == F(3, 2)
    assert internal_potential((0, 0, -1, 0)) == F(-1, 2)
