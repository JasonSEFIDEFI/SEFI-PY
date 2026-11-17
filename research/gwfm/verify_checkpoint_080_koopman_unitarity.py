#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 080
Classical Koopman-von Neumann unitarity.

Repair:
Avoid expand_complex() on transported
real Gaussian exponentials.

No quantum emergence is assumed.
"""

import sympy as sp

x, p, t, s = sp.symbols("x p t s", real=True)
m, a, b = sp.symbols(
    "m a b", positive=True, real=True
)
k = sp.symbols("k", real=True)

I = sp.I
pi = sp.pi

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 080")
print("CLASSICAL KOOPMAN UNITARITY")
print("=" * 64)

H = sp.sqrt(p**2 + m**2)
v = p/H

# 1. Classical measure preservation.
flow = sp.Matrix([x+v*t, p])
J = flow.jacobian([x,p])

assert sp.simplify(J.det()-1) == 0

print("PASS: classical flow preserves phase-space measure")

# 2. Inverse-flow pullback.
# This transports densities forward.
def U(expr, time):
    return expr.subs(
        x, x-v*time, simultaneous=True
    )

f = x**2*p+x*p**2+p**3

assert sp.simplify(
    U(U(f,t),s)-U(f,t+s)
) == 0

assert sp.simplify(U(f,0)-f) == 0

print("PASS: U_s U_t = U_(s+t)")
print("PASS: U_0 is the identity")

# 3. Exact transport equation.
psi0 = sp.exp(-a*x**2/2-b*p**2/2)
psi_t = U(psi0,t)

assert sp.simplify(
    sp.diff(psi_t,t)
    +v*sp.diff(psi_t,x)
) == 0

print("PASS: exact transport evolution")

# 4. Normalized Gaussian.
N = (a*b)**sp.Rational(1,4)/sp.sqrt(pi)

normalized = N*psi0

rho0 = normalized*sp.conjugate(normalized)

expected0 = (
    sp.sqrt(a*b)/pi
    *sp.exp(-a*x**2-b*p**2)
)

assert sp.simplify(rho0-expected0) == 0

norm0 = sp.integrate(
    expected0,
    (x,-sp.oo,sp.oo),
    (p,-sp.oo,sp.oo)
)

assert sp.simplify(norm0-1) == 0

print("PASS: normalized classical L2 state")

# 5. Transported density.
z = x-v*t

transported = U(normalized,t)

rho_t = transported*sp.conjugate(transported)

expected_t = (
    sp.sqrt(a*b)/pi
    *sp.exp(-a*z**2-b*p**2)
)

# Compare without expanding the exponential.
assert sp.simplify(rho_t-expected_t) == 0

# Change variable q=x-v(p)t.
# For each fixed p, dx=dq.
q = sp.symbols("q", real=True)

norm_t = sp.integrate(
    sp.sqrt(a*b)/pi
    *sp.exp(-a*q**2-b*p**2),
    (q,-sp.oo,sp.oo),
    (p,-sp.oo,sp.oo)
)

assert sp.simplify(norm_t-1) == 0

print("PASS: exact L2 norm preservation")

# 6. Classical phase freedom.
psi_A = normalized
psi_B = sp.exp(I*k*x)*normalized

rho_A = psi_A*sp.conjugate(psi_A)
rho_B = psi_B*sp.conjugate(psi_B)

assert sp.simplify(rho_A-rho_B) == 0

# Verify that the functions themselves
# are distinct for a nontrivial phase.
phase_ratio = sp.simplify(
    (psi_B/psi_A).subs({k:1,x:sp.pi})
)

assert phase_ratio == -1

print("PASS: distinct phases encode identical")
print("      classical densities")

# 7. Phase independence after transport.
evolved_A = U(psi_A,t)
evolved_B = U(psi_B,t)

rho_evolved_A = (
    evolved_A*sp.conjugate(evolved_A)
)

rho_evolved_B = (
    evolved_B*sp.conjugate(evolved_B)
)

assert sp.simplify(
    rho_evolved_A-rho_evolved_B
) == 0

print("PASS: classical evolution preserves")
print("      phase-independent densities")

# 8. Classical multiplication operators.
test = x**2+p**2+1

XP = x*(p*test)
PX = p*(x*test)

assert sp.simplify(XP-PX) == 0

print("PASS: classical X and P operators commute")

# 9. Canonical Poisson bracket.
PB_xp = (
    sp.diff(x,x)*sp.diff(p,p)
    -sp.diff(x,p)*sp.diff(p,x)
)

assert PB_xp == 1

print("PASS: {x,p}=1 while [X,P]=0")

# 10. Evolution generator.
generator = -v*sp.diff(psi_t,x)

assert sp.simplify(
    sp.diff(psi_t,t)-generator
) == 0

print("PASS: classical unitary-flow generator")

# 11. Unitarity for general L2 functions
# follows from the measure-preserving,
# invertible flow by change of variables.
# The preceding Gaussian is an exact
# normalized example, not a proof by sampling.

print("\n" + "=" * 64)
print("CHECKPOINT 080: ALL TESTS PASSED")
print("=" * 64)

print("ESTABLISHED MATHEMATICALLY:")
print("- Classical Hamiltonian flow preserves measure.")
print("- Its inverse-flow pullback defines a unitary")
print("  group on complex L2 phase space.")
print("- The classical L2 norm is conserved.")
print("- Different phases can encode the same density.")
print("- Classical multiplication observables commute.")

print("\nCONDITIONAL REPRESENTATION:")
print("- Complex L2 functions were introduced")
print("  as a mathematical representation.")
print("- Squared modulus represents a classical")
print("  phase-space density.")

print("\nNOT ESTABLISHED:")
print("- Physical quantum amplitudes.")
print("- Noncommuting quantum position and momentum.")
print("- A uniquely selected action scale.")
print("- Quantum interference probabilities.")
print("- The quantum Born rule.")
print("- Quantum mechanics emerging from GWFM.")
