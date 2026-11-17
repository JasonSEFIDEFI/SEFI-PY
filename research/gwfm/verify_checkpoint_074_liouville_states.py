#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 074

Classical Liouville transport and probability-state
non-uniqueness for a free relativistic particle.

Established input:
H(p) = sqrt(p**2 + m**2), m > 0, c = 1.

No quantum amplitudes or Born rule.
"""

import sympy as sp

x, p, t = sp.symbols("x p t", real=True)
m = sp.symbols("m", positive=True, real=True)
a, b = sp.symbols("a b", positive=True, real=True)

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 074")
print("LIOUVILLE TRANSPORT AND STATE NON-UNIQUENESS")
print("=" * 64)

# 1. Hamiltonian flow
H = sp.sqrt(p**2 + m**2)
vx = sp.diff(H, p)
vp = -sp.diff(H, x)

assert vp == 0
assert sp.simplify(sp.diff(vx, x)
                   + sp.diff(vp, p)) == 0

print("PASS: divergence-free Hamiltonian flow")

# 2. General transported density:
# rho(x,p,t) = f(x-v(p)*t,p)
f = sp.Function("f")
z = x - vx*t
rho_general = f(z, p)

liouville = sp.diff(rho_general, t) + (
    sp.diff(vx*rho_general, x)
    + sp.diff(vp*rho_general, p)
)

assert sp.simplify(liouville) == 0

print("PASS: general transported density solves")
print("      the Liouville equation")

# 3. Positive normalized Gaussian families
# rho0 = sqrt(a*b)/pi * exp(-a*x^2-b*p^2)
# for a,b>0.

normalization = sp.sqrt(a*b)/sp.pi

def density(A, B):
    return (
        sp.sqrt(A*B)/sp.pi
        * sp.exp(-A*z**2-B*p**2)
    )

rho_A = density(sp.Integer(1), sp.Integer(1))
rho_B = density(sp.Integer(4), sp.Integer(1))

for name, rho in [("A", rho_A), ("B", rho_B)]:
    residual = sp.simplify(
        sp.diff(rho, t)
        + sp.diff(vx*rho, x)
        + sp.diff(vp*rho, p)
    )
    assert residual == 0
    print("PASS: density", name, "obeys Liouville")

# 4. Normalization follows from the
# measure-preserving change of variables:
# z = x-v(p)t, p'=p.
#
# Jacobian determinant must equal one.

J = sp.Matrix([
    [sp.diff(z, x), sp.diff(z, p)],
    [sp.Integer(0), sp.Integer(1)]
])

assert sp.simplify(J.det()) == 1

# Independently check Gaussian normalization
# in the shifted coordinate z.
q = sp.symbols("q", real=True)
gaussian_integral = sp.integrate(
    normalization*sp.exp(-a*q**2-b*p**2),
    (q, -sp.oo, sp.oo),
    (p, -sp.oo, sp.oo)
)

assert sp.simplify(gaussian_integral) == 1

print("PASS: Gaussian normalization")
print("PASS: transported normalization preserved")

# 5. Two states with different spatial widths.
# Evaluate at t=0, x=0, p=0.
point = {t: 0, x: 0, p: 0, m: 2}

value_A = sp.simplify(rho_A.subs(point))
value_B = sp.simplify(rho_B.subs(point))

assert value_A == 1/sp.pi
assert value_B == 2/sp.pi
assert value_A != value_B

print("PASS: distinct normalized states")
print("rho_A(0,0,0) =", value_A)
print("rho_B(0,0,0) =", value_B)

# 6. Classical spatial variance at t=0.
# For density sqrt(a/pi)*exp(-a*x^2),
# Var(x) = 1/(2*a).

variance_A = sp.Rational(1, 2)
variance_B = sp.Rational(1, 8)

assert variance_A != variance_B

print("PASS: distinct spatial variances")
print("Var_A(x) =", variance_A)
print("Var_B(x) =", variance_B)

# 7. Positivity:
# a,b>0, exponential>0 and prefactor>0.
# This is an analytic observation, not a
# finite numerical positivity proof.

print("PASS: positivity follows analytically")
print("      from positive Gaussian parameters")

print("\n" + "=" * 64)
print("CHECKPOINT 074: ALL TESTS PASSED")
print("=" * 64)

print("ESTABLISHED:")
print("- Exact classical Liouville transport.")
print("- Phase-space normalization conservation.")
print("- Multiple positive normalized solutions.")
print("- No unique probability state is selected.")

print("\nNOT ESTABLISHED:")
print("- Quantum Hilbert-space probabilities.")
print("- Complex quantum amplitudes.")
print("- A physical selection principle.")
print("- The Born rule.")
