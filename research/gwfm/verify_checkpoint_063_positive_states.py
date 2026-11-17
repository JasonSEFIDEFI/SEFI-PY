#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 063

Positive states and coherence non-uniqueness.

Purpose:
Demonstrate that two positive normalized states on M_2(C)
can agree on classical history projections but disagree
on an interference-sensitive observable.

Scope:
Finite-dimensional operator algebra.
No physical emergence or Born-rule derivation claimed.
"""

import sympy as sp

half = sp.Rational(1, 2)

I2 = sp.eye(2)

rho_mix = sp.Matrix([
    [half, 0],
    [0, half]
])

rho_coh = sp.Matrix([
    [half, half],
    [half, half]
])

P_A = sp.diag(1, 0)
P_B = sp.diag(0, 1)

H = sp.Matrix([
    [1, 1],
    [1, -1]
]) / sp.sqrt(2)

P_plus = sp.simplify(H.T * P_A * H)


def omega(rho, X):
    return sp.simplify(sp.trace(rho * X))


def check_state(rho):
    assert rho == rho.conjugate().T
    assert sp.trace(rho) == 1

    eigenvalues = list(rho.eigenvals())

    assert all(
        value.is_nonnegative
        for value in eigenvalues
    )

    # For arbitrary X, positivity follows from rho >= 0:
    # Tr(rho X*X) >= 0.
    return eigenvalues


print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 063")
print("POSITIVE STATES AND COHERENCE NON-UNIQUENESS")
print("=" * 64)

print("\nMixture eigenvalues:",
      check_state(rho_mix))

print("Coherent eigenvalues:",
      check_state(rho_coh))

assert omega(rho_mix, I2) == 1
assert omega(rho_coh, I2) == 1

print("PASS: both states normalized and positive")

for P in [P_A, P_B]:
    assert omega(rho_mix, P) == half
    assert omega(rho_coh, P) == half

print("PASS: identical classical history probabilities")

p_mix = omega(rho_mix, P_plus)
p_coh = omega(rho_coh, P_plus)

print("\nIncoherent plus probability:", p_mix)
print("Coherent plus probability:", p_coh)

assert p_mix == half
assert p_coh == 1

print("PASS: different interference outcomes")

# General family of states with identical diagonals.
x, y = sp.symbols("x y", real=True)
z = x + sp.I*y

rho_general = sp.Matrix([
    [half, z],
    [sp.conjugate(z), half]
])

determinant = sp.simplify(
    sp.expand_complex(rho_general.det())
)

print("\nGeneral determinant:", determinant)

assert sp.simplify(
    determinant - (sp.Rational(1, 4) - x**2 - y**2)
) == 0

print("PASS: positivity region x^2 + y^2 <= 1/4")

general_plus = sp.simplify(
    sp.expand_complex(
        omega(rho_general, P_plus)
    )
)

print("General plus probability:", general_plus)

assert sp.simplify(
    general_plus - (half + x)
) == 0

print("PASS: interference depends on off-diagonal coherence")

# An additional measurement accesses imaginary coherence.
P_y = sp.Matrix([
    [half, -sp.I*half],
    [sp.I*half, half]
])

general_y = sp.simplify(
    sp.expand_complex(
        omega(rho_general, P_y)
    )
)

print("Y-basis probability:", general_y)

assert sp.simplify(
    general_y - (half - y)
) == 0

print("PASS: complementary measurement detects imaginary coherence")

print("\n" + "=" * 64)
print("CHECKPOINT 063: ALL TESTS PASSED")
print("=" * 64)

print("CONCLUSION:")
print("Positivity and normalization do not uniquely")
print("determine coherence from classical probabilities.")
print("SCOPE: finite-dimensional M_2(C) state model.")
