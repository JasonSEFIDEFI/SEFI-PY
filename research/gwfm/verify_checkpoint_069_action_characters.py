#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 069

Action additivity versus complex history weights.

Tests:
- Classical action composition
- Continuous exponential character families
- Unit-modulus as an additional constraint
- Non-uniqueness of phase scale

No quantum dynamics or Born rule is derived.
"""

import sympy as sp

S1, S2 = sp.symbols("S1 S2", real=True)
a, b = sp.symbols("a b", real=True)
k = sp.symbols("k", real=True)

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 069")
print("ACTION COMPOSITION AND CHARACTER NON-UNIQUENESS")
print("=" * 64)

# General exponential character.
def F(S):
    return sp.exp((a + sp.I*b)*S)

composition = sp.simplify(
    sp.expand_power_exp(
        F(S1 + S2) / (F(S1)*F(S2))
    )
)

assert composition == 1
assert F(0) == 1

print("PASS: general complex exponential composes")
print("PASS: character normalized at zero action")

# Squared modulus of the general character.
modulus_squared = sp.simplify(
    sp.expand_complex(
        F(S1)*sp.conjugate(F(S1))
    )
)

assert sp.simplify(
    modulus_squared - sp.exp(2*a*S1)
) == 0

print("Squared modulus:", modulus_squared)
print("PASS: nonzero real exponent changes magnitude")

# Explicit counterexamples with identical
# action-composition behavior.
F_unit = sp.exp(sp.I*S1)
F_growth = sp.exp(S1 + sp.I*S1)
F_decay = sp.exp(-S1 + sp.I*S1)

for f in (F_unit, F_growth, F_decay):
    assert sp.simplify(
        f.subs(S1, S1 + S2)
        - f*f.subs(S1, S2)
    ) == 0

print("PASS: unit, growing, and decaying characters compose")

# Distinct unit-modulus characters also compose.
def U(S, scale):
    return sp.exp(sp.I*scale*S)

for scale in (0, 1, 2, sp.Rational(3, 2)):
    assert sp.simplify(
        sp.expand_power_exp(
            U(S1+S2, scale) /
            (U(S1, scale)*U(S2, scale))
        )
    ) == 1

print("PASS: multiple phase scales satisfy composition")

# Explicitly distinguish two phase scales.
test_action = sp.pi/2

phase_one = sp.simplify(U(test_action, 1))
phase_two = sp.simplify(U(test_action, 2))

assert phase_one == sp.I
assert phase_two == -1
assert phase_one != phase_two

print("\nSame geometric action S = pi/2")
print("k = 1 gives:", phase_one)
print("k = 2 gives:", phase_two)
print("PASS: action does not determine phase scale")

# A non-character also exists if we do not
# impose multiplicativity.
G = lambda S: 1 + S**2

failure = sp.expand(
    G(S1+S2) - G(S1)*G(S2)
)

assert failure != 0

print("PASS: action additivity alone does not")
print("      impose multiplicative weights")

print("\n" + "=" * 64)
print("CHECKPOINT 069: ALL TESTS PASSED")
print("=" * 64)

print("CONCLUSION:")
print("Classical geometric action does not uniquely")
print("determine a complex history-weight rule.")
print("Multiplicativity, unit modulus, and the")
print("physical phase scale require justification.")
print("\nNo quantum emergence claim established.")
