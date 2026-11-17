#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 067

Question:
Does exchange symmetry of two histories select
a unique positive normalized state?

Scope:
A conditional two-history matrix representation.
Exchange symmetry is an explicitly imposed
additional condition, not derived from GWFM.
"""

import sympy as sp

half = sp.Rational(1, 2)
x, y = sp.symbols("x y", real=True)

rho = sp.Matrix([
    [half, x + sp.I*y],
    [x - sp.I*y, half]
])

X = sp.Matrix([
    [0, 1],
    [1, 0]
])

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 067")
print("EXCHANGE SYMMETRY AND STATE SELECTION")
print("=" * 64)

# Exchange symmetry
difference = sp.simplify(
    sp.expand_complex(X*rho*X.T - rho)
)

print("\nExchange-symmetry difference:")
sp.pprint(difference)

assert difference == sp.Matrix([
    [0, -2*sp.I*y],
    [2*sp.I*y, 0]
])

print("PASS: exchange symmetry requires y = 0")

# Restrict to exchange-invariant states.
rho_symmetric = rho.subs(y, 0)

assert sp.simplify(
    X*rho_symmetric*X.T - rho_symmetric
) == sp.zeros(2)

eigenvalues = [
    half + x,
    half - x
]

print("\nSymmetric-state eigenvalues:", eigenvalues)
print("Positivity requires -1/2 <= x <= 1/2")

# Three distinct symmetric, positive states.
candidates = {
    "negative coherence": -half,
    "incoherent mixture": sp.Integer(0),
    "positive coherence": half,
}

P_plus = half * sp.Matrix([
    [1, 1],
    [1, 1]
])

def prob(rho_state, P):
    return sp.simplify(sp.trace(rho_state * P))

results = {}

for label, value in candidates.items():
    state = rho_symmetric.subs(x, value)

    assert state == state.conjugate().T
    assert sp.trace(state) == 1
    assert X*state*X.T == state
    assert all(
        e.is_nonnegative
        for e in state.eigenvals()
    )

    results[label] = prob(state, P_plus)

print("\nExchange-invariant states:")
for label, value in results.items():
    print(f"{label:<22} P(+) = {value}")

assert results["negative coherence"] == 0
assert results["incoherent mixture"] == half
assert results["positive coherence"] == 1

print("PASS: three symmetric states give distinct outcomes")

# Even imposing purity leaves two choices.
rho_negative = rho_symmetric.subs(x, -half)
rho_positive = rho_symmetric.subs(x, half)

assert rho_negative**2 == rho_negative
assert rho_positive**2 == rho_positive
assert rho_negative != rho_positive

print("PASS: exchange symmetry plus purity is not unique")

# Symmetry also allows the entire interval.
assert sp.simplify(rho_symmetric.det()) == half**2 - x**2

print("PASS: continuous family of symmetric positive states")

print("\n" + "=" * 64)
print("CHECKPOINT 067: ALL TESTS PASSED")
print("=" * 64)

print("CONCLUSION:")
print("Exchange symmetry alone does not select")
print("a unique coherent positive state.")
print("Even exchange symmetry plus purity leaves")
print("two distinct states in this model.")
print("\nLIMITATIONS:")
print("- The exchange symmetry was imposed.")
print("- The matrix representation was assumed.")
print("- No quantum emergence from geometry is established.")
