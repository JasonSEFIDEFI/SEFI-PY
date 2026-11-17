#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 070
Relativistic worldline action: first and second variations.

Scope:
Classical special relativity, 1+1 dimensions.
No quantum amplitudes or probability assumptions.
"""

import sympy as sp

t = sp.symbols("t", real=True)
T, m = sp.symbols("T m", positive=True, real=True)
u, eps = sp.symbols("u eps", real=True)

# Domain: T > 0, m > 0, |u| < 1.
# Perturbations must remain timelike.

x = u*t + eps*t*(T-t)
velocity = sp.diff(x, t)

L = -m*sp.sqrt(1-velocity**2)

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 070")
print("FIRST AND SECOND VARIATIONS")
print("=" * 64)

# First variation of the Lagrangian.
dL = sp.diff(L, eps).subs(eps, 0)
first = sp.simplify(sp.integrate(dL, (t, 0, T)))

assert first == 0
print("PASS: first variation vanishes")
print("First variation:", first)

# Second variation.
d2L = sp.diff(L, eps, 2).subs(eps, 0)
second = sp.simplify(sp.integrate(d2L, (t, 0, T)))

expected = m*T**3/(3*(1-u**2)**sp.Rational(3, 2))

assert sp.simplify(second-expected) == 0

print("\nSecond variation:")
sp.pprint(second)
print("PASS: second variation positive for |u| < 1")

# Euler-Lagrange equation.
q = sp.Function("q")(t)
qdot = sp.diff(q, t)

L_general = -m*sp.sqrt(1-qdot**2)

momentum = sp.diff(L_general, qdot)
EL = sp.simplify(sp.diff(momentum, t))

expected_EL = (
    m*sp.diff(q, t, 2)
    / (1-sp.diff(q, t)**2)**sp.Rational(3, 2)
)

assert sp.simplify(EL-expected_EL) == 0

print("\nEuler-Lagrange expression:")
sp.pprint(EL)
print("PASS: free stationary paths have zero acceleration")

# Check a straight trajectory.
assert sp.simplify(EL.subs(q, u*t).doit()) == 0
print("PASS: constant-velocity path is stationary")

# Numerical example.
values = {
    T: sp.Integer(2),
    m: sp.Integer(3),
    u: sp.Rational(3, 5)
}

numerical_second = sp.simplify(second.subs(values))

assert numerical_second == sp.Rational(125, 8)

print("\nExact example:")
print("T=2, m=3, u=3/5")
print("Second variation:", numerical_second)
print("PASS: exact numerical verification")

print("\n" + "=" * 64)
print("CHECKPOINT 070: ALL TESTS PASSED")
print("=" * 64)

print("ESTABLISHED:")
print("- Straight timelike trajectories are stationary.")
print("- This perturbation has positive second variation.")
print("- Classical action determines free equations of motion.")

print("\nNOT ESTABLISHED:")
print("- Quantum phase selection.")
print("- Coherent history summation.")
print("- The Born rule.")
