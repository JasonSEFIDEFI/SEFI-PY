"""CP098-T63: Moving adapted-frame dependency audit.

Verify canonical u and q evolution on the regular
constraint surface and determine whether a fixed
adapted frame is preserved.

No gauge generator or physical interpretation assumed.
"""

import sympy as sp

print("=== CP098-T63 ===", flush=True)

v, k = sp.symbols("v k", positive=True)
a, r, s = sp.symbols("a r s", real=True)
C = sp.symbols("C", real=True, nonzero=True)
LT, LN = sp.symbols("lambda_T lambda_N", real=True)

eta = sp.diag(-1, 1, 1, 1)

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])
P = sp.Matrix([0, 0, r, s])

def dot(X, Y):
    return (X.T*eta*Y)[0]

def geom(uu, qq):
    speed = sp.sqrt(-dot(uu, uu))
    T = uu/speed
    w = qq+T*dot(T, qq)
    N = w/sp.sqrt(dot(w,w))
    kap = sp.sqrt(dot(w,w))/speed**2
    return speed, T, N, kap

speed, T, N, kap = geom(u, q)

assert sp.simplify(speed-v) == 0
assert sp.simplify(kap-k) == 0
assert T == sp.Matrix([1, 0, 0, 0])
assert N == sp.Matrix([0, 1, 0, 0])

print("PASS: regular adapted geometry", flush=True)

# Canonical Hamilton equations:
#
# u_dot = dH_T/dp1 = q
#
# q_dot = dH_T/dP
#       = -v^5*k^2*eta*P/(2*C)
#         + LT*T + LN*N

udot = q

qdot = (
    -v**5*k**2*eta*P/(2*C)
    + LT*T
    + LN*N
)

expected_qdot = sp.Matrix([
    LT,
    LN,
    -v**5*k**2*r/(2*C),
    -v**5*k**2*s/(2*C)
])

assert qdot == expected_qdot

print("PASS: canonical external evolution",
      flush=True)

print("u_dot =", udot, flush=True)
print("q_dot =", qdot, flush=True)

# A fixed adapted frame requires the spatial
# components of u to remain zero.
#
# But u_dot[1] = v^2*k is nonzero throughout
# the regular sector v>0, k>0.

assert sp.factor(udot[1]-v**2*k) == 0

print("PASS: fixed adapted frame is not preserved",
      flush=True)

# Compute the actual tangent-vector evolution.
#
# T=u/v, so
#
# T_dot=(u_dot+T*(T.u_dot))/v.
#
# The projection removes the component parallel
# to T and automatically accounts for v_dot.

Tdot = sp.simplify(
    (udot+T*dot(T,udot))/v
)

expected_Tdot = sp.Matrix([0,v*k,0,0])

assert Tdot == expected_Tdot

print("T_dot =", Tdot, flush=True)
print("PASS: tangent-frame transport", flush=True)

# Verify normalization is preserved to
# first order: T.T=-1 implies T.T_dot=0.

assert sp.simplify(dot(T,Tdot)) == 0

print("PASS: tangent normalization preserved",
      flush=True)

# A phase-space function evaluated in an
# adapted frame generally requires both
# intrinsic derivative terms and derivatives
# arising from frame transport.
#
# T60 computed E at an adapted point.
# Its printed expression is not automatically
# a global coordinate expression.

print("ADAPTED FRAME: EVOLVES NONTRIVIALLY")
print("T_DOT = v*k*N: VERIFIED")
print("FIXED-FRAME DIFFERENTIATION: INSUFFICIENT")
print("FULL CANONICAL E EVOLUTION: OPEN")
print("FRAME-TRANSPORT CONTRIBUTIONS: OPEN")
print("GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
