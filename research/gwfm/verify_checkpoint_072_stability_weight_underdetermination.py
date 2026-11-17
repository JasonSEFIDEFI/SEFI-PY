#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 072

Does classical variational stability uniquely determine
a complex history-weight prescription?

Established input:
S[x] = -m integral sqrt(1-xdot**2) dt

Conditional constructions:
W_k(S) = exp(i*k*S)
W_d(S) = exp(-S)

These are mathematical counterexamples to uniqueness,
not derived quantum amplitudes.
"""

import sympy as sp

t = sp.symbols("t", real=True)
T, m = sp.symbols("T m", positive=True, real=True)
u, eps = sp.symbols("u eps", real=True)
S1, S2 = sp.symbols("S1 S2", real=True)

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 072")
print("CLASSICAL STABILITY VS HISTORY-WEIGHT SELECTION")
print("=" * 64)

# --------------------------------------------------
# 1. Classical geometric action
# --------------------------------------------------

eta = t*(T-t)
velocity = u + eps*sp.diff(eta, t)

L = -m*sp.sqrt(1-velocity**2)

first = sp.simplify(
    sp.integrate(
        sp.diff(L, eps).subs(eps, 0),
        (t, 0, T)
    )
)

second = sp.simplify(
    sp.integrate(
        sp.diff(L, eps, 2).subs(eps, 0),
        (t, 0, T)
    )
)

expected = (
    m*T**3 /
    (3*(1-u**2)**sp.Rational(3, 2))
)

assert first == 0
assert sp.simplify(second-expected) == 0

print("PASS: classical stationary trajectory")
print("PASS: positive second-variation formula")

# --------------------------------------------------
# 2. Alternative history weights
# --------------------------------------------------

def W1(S):
    return sp.exp(sp.I*S)

def W2(S):
    return sp.exp(2*sp.I*S)

def W3(S):
    return sp.exp(-S)

weights = {
    "unit phase k=1": W1,
    "unit phase k=2": W2,
    "real exponential": W3,
}

for name, W in weights.items():
    composition = sp.simplify(
        sp.expand_power_exp(
            W(S1+S2)/(W(S1)*W(S2))
        )
    )
    assert composition == 1
    assert W(0) == 1
    print("PASS: composition:", name)

# --------------------------------------------------
# 3. Same geometric action, different weights
# --------------------------------------------------

action_value = sp.pi/2

results = {
    name: sp.simplify(W(action_value))
    for name, W in weights.items()
}

print("\nSame action S=pi/2:")

for name, value in results.items():
    print(name, "->", value)

assert results["unit phase k=1"] == sp.I
assert results["unit phase k=2"] == -1
assert results["real exponential"] == sp.exp(-sp.pi/2)

assert len(set(results.values())) == 3

print("PASS: same action admits distinct weights")

# --------------------------------------------------
# 4. Local stability does not distinguish weights
# --------------------------------------------------

# The classical Hessian is independent of the
# subsequently chosen weighting prescription.

hessian = second

for name in weights:
    assert sp.simplify(hessian-expected) == 0

print("PASS: classical Hessian unchanged by weight choice")

# --------------------------------------------------
# 5. Unit modulus still leaves ambiguity
# --------------------------------------------------

for W in (W1, W2):
    modulus_squared = sp.simplify(
        sp.expand_complex(
            W(S1)*sp.conjugate(W(S1))
        )
    )
    assert modulus_squared == 1

assert sp.simplify(
    W1(sp.pi/2)-W2(sp.pi/2)
) != 0

print("PASS: unit modulus does not select phase scale")

print("\n" + "=" * 64)
print("CHECKPOINT 072: ALL TESTS PASSED")
print("=" * 64)

print("RESULT:")
print("Classical action composition and local")
print("variational stability do not uniquely select")
print("a complex history-weight prescription.")

print("\nSTOP-GATE:")
print("An additional justified selection principle")
print("is required before a physical quantum")
print("history-weight law can be claimed.")
