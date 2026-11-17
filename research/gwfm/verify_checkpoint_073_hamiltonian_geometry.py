#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 073 (corrected)

Relativistic Hamiltonian geometry.

Physical assumptions:
m > 0, |v| < 1, c = 1.

Classical results only.
"""

import sympy as sp

v, p, t = sp.symbols("v p t", real=True)
m = sp.symbols("m", positive=True, real=True)
x, p0, x0 = sp.symbols("x p0 x0", real=True)

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 073")
print("RELATIVISTIC HAMILTONIAN GEOMETRY")
print("=" * 64)

# 1. Relativistic Lagrangian
L = -m*sp.sqrt(1-v**2)

canonical_p = sp.diff(L, v)

assert sp.simplify(
    canonical_p-m*v/sp.sqrt(1-v**2)
) == 0

print("PASS: canonical momentum derived")

# 2. Invert momentum-velocity relation
velocity_of_p = p/sp.sqrt(p**2+m**2)

assert sp.simplify(
    canonical_p.subs(v, velocity_of_p)-p
) == 0

print("PASS: momentum-velocity inversion")

# 3. Legendre transform
H_legendre = sp.simplify(
    p*velocity_of_p-L.subs(v, velocity_of_p)
)

H = sp.sqrt(p**2+m**2)

assert sp.simplify(H_legendre-H) == 0

print("PASS: Hamiltonian derived")
print("H =", H)

# 4. Hamilton's equations
xdot = sp.diff(H, p)
pdot = -sp.diff(H, x)

assert sp.simplify(xdot-velocity_of_p) == 0
assert pdot == 0

print("PASS: Hamilton's equations")

# 5. Liouville divergence
divergence = sp.simplify(
    sp.diff(xdot, x)+sp.diff(pdot, p)
)

assert divergence == 0

print("PASS: phase-space divergence vanishes")

# 6. Exact Hamiltonian flow
x_t = x0+t*p0/sp.sqrt(p0**2+m**2)
p_t = p0

J = sp.Matrix([
    [sp.diff(x_t, x0), sp.diff(x_t, p0)],
    [sp.diff(p_t, x0), sp.diff(p_t, p0)]
])

assert sp.simplify(J.det()) == 1

print("PASS: phase-space area preserved")

# 7. Canonical Poisson bracket
def poisson(f, g):
    return sp.simplify(
        sp.diff(f, x0)*sp.diff(g, p0)
        - sp.diff(f, p0)*sp.diff(g, x0)
    )

assert poisson(x_t, p_t) == 1

print("PASS: canonical Poisson bracket preserved")

# 8. Symplectic form
Omega = sp.Matrix([[0, 1], [-1, 0]])

assert sp.simplify(
    J.T*Omega*J-Omega
) == sp.zeros(2)

print("PASS: symplectic form preserved")

# 9. Exact numerical example
values = {
    m: sp.Integer(2),
    p0: sp.Rational(3, 2),
    x0: sp.Integer(0),
    t: sp.Integer(5)
}

x_final = sp.simplify(x_t.subs(values))
p_final = sp.simplify(p_t.subs(values))

assert x_final == 3
assert p_final == sp.Rational(3, 2)

print("PASS: exact numerical flow")
print("x(5) =", x_final)
print("p(5) =", p_final)

print("\n" + "=" * 64)
print("CHECKPOINT 073: ALL TESTS PASSED")
print("=" * 64)

print("ESTABLISHED:")
print("- Relativistic Hamiltonian structure.")
print("- Canonical symplectic evolution.")
print("- Phase-space area preservation.")
print("- Canonical Poisson-bracket preservation.")

print("\nNOT ESTABLISHED:")
print("- Quantum Hilbert-space structure.")
print("- Quantum unitarity.")
print("- Quantum amplitudes or Born probabilities.")
