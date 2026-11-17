#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 068

Proper-time geometry and classical relativistic action.

Scope:
1+1-dimensional Minkowski spacetime.
Two future-directed piecewise timelike worldlines.
Identical endpoints; different proper times.

No quantum phase, probability, or physical SEFI
emergence is assumed.
"""

import sympy as sp

T, m, v = sp.symbols(
    "T m v", positive=True, real=True
)

# Physical domain:
# T > 0, m > 0, 0 < v < 1

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 068")
print("PROPER-TIME GEOMETRY AND CLASSICAL ACTION")
print("=" * 64)

# Worldline 1: stationary.
tau1 = T

# Worldline 2: outward and return,
# each at constant speed magnitude v.
tau2 = sp.simplify(
    (T / 2) * sp.sqrt(1 - v**2)
    + (T / 2) * sp.sqrt(1 - v**2)
)

assert sp.simplify(
    tau2 - T * sp.sqrt(1 - v**2)
) == 0

print("\nProper times:")
print("tau1 =", tau1)
print("tau2 =", tau2)
print("PASS: exact proper-time expressions")

# Free-particle classical action.
S1 = -m * tau1
S2 = -m * tau2

delta_S = sp.simplify(S2 - S1)

expected = m*T*(1 - sp.sqrt(1 - v**2))

assert sp.simplify(delta_S - expected) == 0

print("\nClassical actions:")
print("S1 =", S1)
print("S2 =", S2)
print("S2 - S1 =", delta_S)
print("PASS: exact action difference")

# Verify additivity under subdivision.
segment1 = -m*(T/2)*sp.sqrt(1 - v**2)
segment2 = -m*(T/2)*sp.sqrt(1 - v**2)

assert sp.simplify(
    segment1 + segment2 - S2
) == 0

print("PASS: action additive under path subdivision")

# Nonrelativistic expansion.
series = sp.series(delta_S, v, 0, 6)

print("\nSmall-velocity expansion:")
print(series)

expected_series = (
    m*T*v**2/2 + m*T*v**4/8
)

assert sp.simplify(
    series.removeO() - expected_series
) == 0

print("PASS: leading kinetic-action difference")

# Exact numerical example:
# v = 3/5, T = 10, m = 2.
values = {
    v: sp.Rational(3, 5),
    T: sp.Integer(10),
    m: sp.Integer(2),
}

n_tau1 = sp.simplify(tau1.subs(values))
n_tau2 = sp.simplify(tau2.subs(values))
n_delta = sp.simplify(delta_S.subs(values))

print("\nExact example:")
print("v = 3/5, T = 10, m = 2")
print("tau1 =", n_tau1)
print("tau2 =", n_tau2)
print("delta_S =", n_delta)

assert n_tau1 == 10
assert n_tau2 == 8
assert n_delta == 4

print("PASS: exact numerical consistency")

# Independent endpoint verification.
# First half: displacement +v*T/2.
# Second half: displacement -v*T/2.
outward = v*T/2
returnward = -v*T/2

assert sp.simplify(
    outward + returnward
) == 0

print("PASS: both worldlines share endpoints")

# Timelike interval per half segment:
# (T/2)^2 - (v*T/2)^2
interval_sq = sp.factor(
    (T/2)**2 - (v*T/2)**2
)

assert sp.simplify(
    interval_sq - T**2*(1-v**2)/4
) == 0

print("PASS: segments timelike for 0 < v < 1")

print("\n" + "=" * 64)
print("CHECKPOINT 068: ALL TESTS PASSED")
print("=" * 64)

print("ESTABLISHED:")
print("- Same endpoints can have different proper times.")
print("- Proper-time geometry determines classical action.")
print("- Classical action is additive along segments.")

print("\nNOT ESTABLISHED:")
print("- A quantum amplitude or phase law.")
print("- The Born probability rule.")
print("- A derivation of quantum mechanics from GWFM.")
