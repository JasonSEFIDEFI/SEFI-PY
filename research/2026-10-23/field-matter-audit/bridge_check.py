"""Limited symmetry bridge audit, not a derivation of a physical model."""
from fractions import Fraction as F
import itertools
import json
from pathlib import Path


def invariant_monomials(max_degree=4):
    # Powers of Phi, Phi*, Sigma, Sigma*. Independent U(1) rotations.
    return [powers for powers in itertools.product(range(max_degree+1), repeat=4)
            if sum(powers) <= max_degree and powers[0] == powers[1] and powers[2] == powers[3]]


def portal_identity(a, b, c, g, v2):
    """General portal potential with canonical real-component normalization.

    U = a rho^2/4+b(eta-v2)^2/4+c rho/2+g rho eta/2.
    U-m^2 rho/2 = b(eta-v2+g rho/b)^2/4+(ab-g^2)rho^2/(4b).
    """
    a, b, c, g, v2 = map(F, (a, b, c, g, v2))
    mass2 = c+g*v2
    lhs = {(2, 0): a/4, (0, 2): b/4, (0, 1): -b*v2/2,
           (0, 0): b*v2*v2/4, (1, 0): (c-mass2)/2, (1, 1): g/2}
    rhs = {(2, 0): g*g/(4*b)+(a*b-g*g)/(4*b), (0, 2): b/4,
           (0, 1): -b*v2/2, (0, 0): b*v2*v2/4,
           (1, 0): -g*v2/2, (1, 1): g/2}
    return {'identity_exact': lhs == rhs, 'm_squared': str(mass2),
            'quartic_discriminant_ab_minus_g2': str(a*b-g*g)}


def internal_potential(I):
    # Allowed constant-coupling example from the foundations worldline action:
    # all stiffnesses=1, mu1=1, kappa=1, torsion=0.
    return sum(F(x)**2 for x in I)/2+F(I[2])


def findings():
    return {'scope': 'Polynomial potential operators only, not all EFT derivative operators',
            'invariant_powers_phi_phistar_sigma_sigmastar': invariant_monomials(),
            'lambda12': portal_identity(1, 12, 1, 4, F(1, 4)),
            'lambda25': portal_identity(1, 25, 1, 4, F(3, 25)),
            'naive_internal_pairing_counterexample': {
                'pairing': 'Phi=I_O+i I_A; Sigma=I_S+i I_W',
                'state': [0, 0, 1, 0], 'rotated_state': [0, 0, -1, 0],
                'rotation': 'Sigma -> -Sigma',
                'potential_before': str(internal_potential((0, 0, 1, 0))),
                'potential_after': str(internal_potential((0, 0, -1, 0))),
                'interpretation': 'This allowed worldline coupling lacks the bulk U(1)xU(1) symmetry under naive pairing. Not a no-go theorem for every possible bridge.'}}


if __name__ == '__main__':
    path = Path(__file__).with_name('bridge-results.json')
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(findings(), indent=2)+'\n')
