"""CP098-T68: External scalar transport audit.

Verify total derivatives of geometric contractions
on the regular constraint surface.

This is a dependency test, not the full E bracket.
"""

import sympy as sp

print("=== CP098-T68 ===", flush=True)

v, k = sp.symbols("v k", positive=True)
a, r, s, x, y, n, b, c = sp.symbols(
    "a r s x y n b c", real=True
)
A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)
LT, LN = sp.symbols(
    "lambda_T lambda_N", real=True
)

eta = sp.diag(-1, 1, 1, 1)

def dot(X, Y):
    return (X.T * eta * Y)[0]

def pair(p, X):
    return (p.T * X)[0]

def eq(X, Y):
    return sp.cancel(sp.together(X-Y)) == 0

def derivative_pair(p, pdot, X, Xdot):
    return sp.factor(
        pair(pdot, X) + pair(p, Xdot)
    )

T = sp.Matrix([1, 0, 0, 0])
N = sp.Matrix([0, 1, 0, 0])

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])
P = sp.Matrix([0, 0, r, s])
p1 = sp.Matrix([0, -2*A*k/v, x, y])

R = r**2+s**2
S = r*x+s*y

p0 = sp.Matrix([
    A*k**2+v**4*k**2*R/(4*C)-B,
    n, b, c
])

udot = q

qdot = sp.Matrix([
    LT,
    LN,
    -v**5*k**2*r/(2*C),
    -v**5*k**2*s/(2*C)
])

Tdot = sp.Matrix([0, v*k, 0, 0])

Ndot = sp.Matrix([
    v*k,
    0,
    -v**3*k*r/(2*C),
    -v**3*k*s/(2*C)
])

Pdot = sp.Matrix([
    0,
    v**3*k*R/(2*C),
    -x-LN*r/(v**2*k),
    -y-LN*s/(v**2*k)
])

p1dot = sp.Matrix([
    2*A*k**2,
    2*A*a*k/v**2-n-a*k*v**2*R/(2*C),
    a*LN*r/(k*v**3)-LT*r/v-b,
    a*LN*s/(k*v**3)-LT*s/v-c
])

p0dot = sp.zeros(4, 1)

# Check basic adapted constraints.
assert eq(pair(P,T), 0)
assert eq(pair(P,N), 0)
assert eq(pair(p1,T), 0)
assert eq(pair(p1,N), -2*A*k/v)

print("PASS: regular adapted constraints", flush=True)

# Primary normal constraint: two contributions
# must cancel.
phiNdot = derivative_pair(P, Pdot, N, Ndot)

assert eq(phiNdot, 0)

print("PASS: phi_N transport cancellation",
      flush=True)

# Tangential secondary constraint.
chiTdot = derivative_pair(p1, p1dot, T, Tdot)

assert eq(chiTdot, 0)

print("PASS: chi_T transport cancellation",
      flush=True)

# Tangential primary constraint:
# P_dot.T + P.T_dot.
phiTdot = derivative_pair(P, Pdot, T, Tdot)

assert eq(phiTdot, 0)

print("PASS: phi_T transport", flush=True)

# Geometric acceleration component
# a = -T.q.
adot = -(
    dot(Tdot, q) + dot(T, qdot)
)

expected_adot = LT-v**3*k**2

assert eq(adot, expected_adot)

print("PASS: a transport", flush=True)
print("a_dot =", sp.factor(adot), flush=True)

# n = p0.N. Although p0_dot=0,
# the moving normal contributes.
ndot = derivative_pair(p0, p0dot, N, Ndot)

expected_ndot = (
    v*k*p0[0]
    - v**3*k*(b*r+c*s)/(2*C)
)

assert eq(ndot, expected_ndot)

print("PASS: p0.N transport", flush=True)
print("n_dot =", sp.factor(ndot), flush=True)

# The Minkowski momentum norm is a scalar.
# Its derivative is 2 P.P_dot.
Rdot = sp.factor(2*dot(P, eta*Pdot))

# Here eta*Pdot converts covector derivative
# to its corresponding contravariant vector.
expected_Rdot = sp.factor(
    -2*S-2*LN*R/(v**2*k)
)

assert eq(Rdot, expected_Rdot)

print("PASS: P_squared transport", flush=True)
print("P_squared_dot =", Rdot, flush=True)

# The full contraction P_mu p1^mu is
# invariantly defined. It must be
# differentiated with the Minkowski metric.
Sdot = sp.factor(
    pair(Pdot, eta*p1)
    + pair(P, eta*p1dot)
)

expected_Sdot = sp.factor(
    -x**2-y**2-r*b-s*c
    -LN*S/(v**2*k)
    -LT*R/v
    +a*LN*R/(k*v**3)
    -A*k**2*R*v**2/C
)

assert eq(Sdot, expected_Sdot)

print("PASS: P.p1 transport", flush=True)
print("P.p1_dot =", Sdot, flush=True)

# The transverse components of P and p1
# define S=r*x+s*y on the adapted slice.
# Their coordinate derivative is not
# generally the same as the invariant
# contraction derivative when the frame moves.
#
# Verify the distinction explicitly.

coordinate_Sdot = sp.factor(
    Pdot[2]*x + Pdot[3]*y
    +r*p1dot[2]+s*p1dot[3]
)

difference = sp.factor(
    Sdot-coordinate_Sdot
)

expected_difference = -A*k**2*R*v**2/C

assert eq(difference, expected_difference)

print(
    "Invariant minus coordinate S_dot =",
    difference,
    flush=True
)

print(
    "PASS: invariant/coordinate correction",
    flush=True
)

# No equality is assumed between these two
# derivatives; this guards against treating
# adapted coordinates as fixed global scalars.

print("PASS: invariant/coordinate distinction retained",
      flush=True)

print("EXTERNAL SCALAR TRANSPORT: VERIFIED")
print("FULL E EVOLUTION: OPEN")
print("JACOBI CROSS-CHECK: OPEN")
print("GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
