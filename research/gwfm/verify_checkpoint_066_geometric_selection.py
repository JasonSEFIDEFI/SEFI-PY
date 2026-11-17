#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 066

Question:
Does fixed geometric action, even together with a
specified action-phase transformation, uniquely
select a coherent positive state?

This is a counterexample test in a conditional
two-history quantum comparison model.

No derivation of quantum mechanics is claimed.
"""

import sympy as sp

half = sp.Rational(1, 2)

# Fixed action data for two distinct histories.
S1 = sp.Integer(0)
S2 = sp.pi

# A chosen phase scale. This is additional structure.
k = sp.Integer(1)

D = sp.diag(
    sp.exp(sp.I * k * S1),
    sp.exp(sp.I * k * S2)
)

D = sp.simplify(D)

# Same classical history probabilities.
rho_mix = half * sp.eye(2)

rho_coh = sp.Matrix([
    [half, half],
    [half, half]
])

def transform(rho):
    return sp.simplify(D * rho * D.conjugate().T)

def probability(rho, projector):
    return sp.simplify(sp.trace(rho * projector))

P_A = sp.diag(1, 0)
P_B = sp.diag(0, 1)

P_plus = half * sp.Matrix([
    [1, 1],
    [1, 1]
])

def is_state(rho):
    assert rho == rho.conjugate().T
    assert sp.trace(rho) == 1
    assert all(
        e.is_nonnegative
        for e in rho.eigenvals()
    )

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 066")
print("GEOMETRIC ACTION AND STATE NON-UNIQUENESS")
print("=" * 64)

print("\nFixed actions:", S1, S2)
print("Chosen phase scale:", k)
print("Phase transformation:")
sp.pprint(D)

for rho in (rho_mix, rho_coh):
    is_state(rho)
    is_state(transform(rho))

print("\nPASS: both states remain positive and normalized")

for rho in (rho_mix, rho_coh):
    assert probability(rho, P_A) == half
    assert probability(rho, P_B) == half
    assert probability(transform(rho), P_A) == half
    assert probability(transform(rho), P_B) == half

print("PASS: identical classical history probabilities")

mix_after = transform(rho_mix)
coh_after = transform(rho_coh)

p_mix = probability(mix_after, P_plus)
p_coh = probability(coh_after, P_plus)

print("\nRecombined plus probabilities:")
print("Mixture:", p_mix)
print("Coherent:", p_coh)

assert p_mix == half
assert p_coh == 0

print("PASS: same actions permit different outcomes")

# The phase transformation preserves coherence
# magnitude. It does not create or select it.
assert mix_after[0, 1] == 0
assert sp.simplify(
    coh_after[0, 1] *
    sp.conjugate(coh_after[0, 1])
) == sp.Rational(1, 4)

print("PASS: action-phase map does not select coherence")

# Continuum of admissible initial coherence values.
r = sp.symbols("r", real=True)

rho_r = sp.Matrix([
    [half, r],
    [r, half]
])

assert sp.simplify(rho_r.det()) == half**2 - r**2

print("\nGeneral family:")
print("rho(r) eigenvalues: 1/2 + r, 1/2 - r")
print("Positivity: -1/2 <= r <= 1/2")

assert sp.simplify(
    probability(transform(rho_r), P_plus)
    - (half - r)
) == 0

print("Recombined plus probability: 1/2 - r")

print("\n" + "=" * 64)
print("CHECKPOINT 066: ALL TESTS PASSED")
print("=" * 64)
print("Conditional conclusion:")
print("Fixed actions and the chosen phase map do not")
print("uniquely determine the initial positive state.")
print("This does not rule out stronger geometric laws.")
