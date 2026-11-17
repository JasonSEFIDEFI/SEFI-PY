#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 079

Topology of the established 1+1-dimensional
free relativistic particle phase space.

Established model:
M = R^2 with coordinates (x,p)
omega = dx wedge dp
H = sqrt(p^2 + m^2), m > 0

Distinguish:
- Exact symplectic geometry on R^2
- Conditional prequantization on T^2
- Physical quantum emergence (not established)
"""

import sympy as sp

x, p, t = sp.symbols("x p t", real=True)
m, R, B = sp.symbols(
    "m R B", positive=True, real=True
)
lam = sp.symbols("lambda", real=True)
s = sp.symbols("s", real=True)

Omega = sp.Matrix([[0, 1], [-1, 0]])

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 079")
print("TOPOLOGY OF RELATIVISTIC PHASE SPACE")
print("=" * 64)

# 1. The global one-form theta = x dp
# satisfies d(theta) = dx wedge dp.
#
# Represent a one-form by its dx, dp coefficients.
theta_x = sp.Integer(0)
theta_p = x

dtheta = sp.diff(theta_p, x) - sp.diff(theta_x, p)

assert dtheta == 1
print("PASS: omega = d(x dp) globally on R^2")

# 2. The alternative canonical one-form p dx
# has derivative -omega.
alpha_x = p
alpha_p = sp.Integer(0)

dalpha = sp.diff(alpha_p, x) - sp.diff(alpha_x, p)

assert dalpha == -1

# Their difference is exact:
# x dp + p dx = d(xp).
assert sp.diff(x*p, x) == p
assert sp.diff(x*p, p) == x

print("PASS: canonical one-form sign conventions consistent")

# 3. Relativistic Hamiltonian vector field.
H = sp.sqrt(p**2 + m**2)

velocity = sp.diff(H, p)
acceleration = -sp.diff(H, x)

assert sp.simplify(
    velocity - p/H
) == 0

assert acceleration == 0

print("PASS: established relativistic Hamiltonian flow")

# 4. Time-t flow:
# Phi_t(x,p) = (x + v(p)*t, p).
flow_x = x + velocity*t
flow_p = p

J = sp.Matrix([
    [sp.diff(flow_x,x), sp.diff(flow_x,p)],
    [sp.diff(flow_p,x), sp.diff(flow_p,p)]
])

assert sp.simplify(J.det()-1) == 0
assert sp.simplify(
    J.T*Omega*J-Omega
) == sp.zeros(2)

print("PASS: free flow preserves the symplectic form")

# 5. Verify the pullback of theta=x dp.
# Phi_t^*theta = (x+v*t) dp.
# Difference from theta is t*v(p) dp.
#
# It is exact, with primitive t*H(p).
assert sp.simplify(
    (flow_x-x)-sp.diff(t*H,p)
) == 0

print("PASS: flow changes theta by an exact one-form")

# 6. A square loop in R^2 is contractible.
# Its p dx integral equals minus its area.

C = [
    (sp.Integer(0),sp.Integer(0)),
    (R,sp.Integer(0)),
    (R,R),
    (sp.Integer(0),R),
    (sp.Integer(0),sp.Integer(0))
]

def loop_integral(vertices):
    total = sp.Integer(0)
    for (x1,p1),(x2,p2) in zip(
        vertices[:-1],vertices[1:]
    ):
        total += (p1+p2)*(x2-x1)/2
    return sp.simplify(total)

action = loop_integral(C)

assert sp.simplify(action+R**2) == 0

print("PASS: contractible loop obeys Stokes' theorem")

# 7. Explicit contraction:
# C_s = s*C, 0 <= s <= 1.
# Its enclosed area scales as s^2.
contracted = [(s*a,s*b) for a,b in C]

contracted_action = loop_integral(contracted)

assert sp.simplify(
    contracted_action-s**2*action
) == 0

assert sp.simplify(
    contracted_action.subs(s,0)
) == 0

assert sp.simplify(
    contracted_action.subs(s,1)-action
) == 0

print("PASS: loop contracts continuously to a point")

# 8. Explicit global primitive for omega
# demonstrates exactness on the plane.
#
# For an exact two-form, its integral over
# any compact, closed, oriented 2-cycle is zero
# by Stokes' theorem.
#
# R^2 has H_2(R^2,Z)=0; this topological
# statement is a standard theorem, not
# established by a finite symbolic test.

assert dtheta == 1
print("PASS: global exactness verified symbolically")
print("NOTE: H_2(R^2,Z)=0 is an external theorem")

# 9. On the torus, omega=B dx wedge dy
# has total integral B.
#
# Under an additional U(1) prequantization
# requirement, lambda*B/(2*pi) is integral.
n = sp.symbols("n", integer=True)

allowed_lambda = 2*sp.pi*n/B

assert sp.simplify(
    allowed_lambda*B/(2*sp.pi)-n
) == 0

print("PASS: torus integrality remains conditional")

# 10. On R^2, the curvature prescription
# i*lambda*omega can be implemented using
# the globally defined connection one-form
# i*lambda*x dp for every real lambda.
#
# This is a mathematical construction,
# not a physical selection rule.
connection_coefficient = lam*x
curvature = sp.diff(connection_coefficient,x)

assert sp.simplify(curvature-lam) == 0

print("PASS: global plane connection exists")
print("      for arbitrary real lambda")

# 11. Two different positive scales both work.
lambda_a = sp.Integer(1)
lambda_b = sp.Integer(2)

assert lambda_a != lambda_b
assert sp.diff(lambda_a*x,x) == lambda_a
assert sp.diff(lambda_b*x,x) == lambda_b

print("PASS: plane geometry does not select lambda")

print("\n" + "=" * 64)
print("CHECKPOINT 079: ALL TESTS PASSED")
print("=" * 64)

print("ESTABLISHED IN THE CLASSICAL MODEL:")
print("- The phase space is R^2.")
print("- Its symplectic form has a global primitive.")
print("- The free relativistic flow is symplectic.")
print("- Contractible loops carry geometric action.")
print("- The tested loop contracts to a point.")

print("\nSTANDARD TOPOLOGICAL INPUT:")
print("- R^2 has no nontrivial compact 2-cycles.")
print("- The torus has a nontrivial 2-cycle.")

print("\nCONDITIONAL CONSTRUCTION:")
print("- A global U(1) connection with curvature")
print("  i*lambda*omega exists on R^2 for")
print("  arbitrary real lambda.")

print("\nNOT ESTABLISHED:")
print("- A toroidal GWFM phase space.")
print("- A physically required U(1) connection.")
print("- A unique quantum action scale.")
print("- Emergent quantum mechanics.")
print("- The Born rule.")
