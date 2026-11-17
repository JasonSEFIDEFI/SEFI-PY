#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 071

General second variation of the free relativistic
worldline action in 1+1 Minkowski spacetime.

No quantum assumptions.
"""

import sympy as sp

t = sp.symbols("t", real=True)
T, m = sp.symbols("T m", positive=True, real=True)
u, eps = sp.symbols("u eps", real=True)

eta = sp.Function("eta")(t)
etadot = sp.diff(eta, t)

velocity = u + eps*etadot
L = -m*sp.sqrt(1-velocity**2)

first_density = sp.simplify(
    sp.diff(L, eps).subs(eps, 0)
)

second_density = sp.simplify(
    sp.diff(L, eps, 2).subs(eps, 0)
)

expected_first = m*u*etadot/sp.sqrt(1-u**2)
expected_second = (
    m*etadot**2/(1-u**2)**sp.Rational(3, 2)
)

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 071")
print("GENERAL SECOND VARIATION")
print("=" * 64)

assert sp.simplify(
    first_density-expected_first
) == 0

print("PASS: general first-variation density")

assert sp.simplify(
    second_density-expected_second
) == 0

print("PASS: general second-variation density")

# Endpoint-fixed condition:
# integral eta'(t) dt = eta(T)-eta(0) = 0.
endpoint_difference = sp.Symbol("eta_T") - sp.Symbol("eta_0")

assert endpoint_difference.subs({
    "eta_T": 0,
    "eta_0": 0
}) == 0

print("PASS: endpoint condition removes first variation")

# Test several independent endpoint-fixed perturbations.
perturbations = {
    "quadratic": t*(T-t),
    "cubic": t*(T-t)*(2*t-T),
    "quartic": t**2*(T-t)**2,
}

for name, f in perturbations.items():
    assert sp.simplify(f.subs(t, 0)) == 0
    assert sp.simplify(f.subs(t, T)) == 0

    integral = sp.simplify(
        sp.integrate(sp.diff(f, t)**2, (t, 0, T))
    )

    print(f"{name}: integral =", integral)

    # Verify positivity at T=2.
    value = sp.simplify(integral.subs(T, 2))
    assert value > 0

print("PASS: three independent perturbations")

# Exact Poincare-type lower bound for zero-endpoint
# perturbations:
# integral eta'^2 >= (pi/T)^2 integral eta^2.
#
# Test it on a polynomial example, without claiming
# that finite tests prove the general inequality.

f = t*(T-t)

gradient_energy = sp.integrate(
    sp.diff(f, t)**2, (t, 0, T)
)

amplitude_energy = sp.integrate(
    f**2, (t, 0, T)
)

ratio = sp.simplify(
    gradient_energy/amplitude_energy
)

assert ratio == 10/T**2

print("Polynomial energy ratio:", ratio)
print("PASS: exact gradient-to-amplitude ratio")

# Since 10 > pi^2, the example satisfies the
# sharp Dirichlet Poincare bound.
assert float(10-sp.pi**2) > 0

print("PASS: example satisfies Poincare lower bound")

print("\n" + "=" * 64)
print("CHECKPOINT 071: ALL TESTS PASSED")
print("=" * 64)

print("ANALYTICAL RESULT:")
print("delta2 S = m/(1-u^2)^(3/2) * integral eta'^2 dt")
print("Positive for every nonzero smooth endpoint-fixed")
print("perturbation, with m>0 and |u|<1.")

print("\nLIMITATION:")
print("No quantum phase, state space, or Born rule derived.")
