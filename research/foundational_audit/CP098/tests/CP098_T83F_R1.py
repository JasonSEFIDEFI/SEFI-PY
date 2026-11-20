"""CP098-T83F-R1: Independent canonical n=p0.N transport.

Derive Ndot from canonical u/q velocities.
Apply regular constraints only after differentiation.
No T74 transport equations imported.
"""

import sympy as sp

print("=== CP098-T83F-R1 ===", flush=True)

v, k = sp.symbols("v k", positive=True)
a, LT, LN = sp.symbols("a LT LN", real=True)
A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)

eps = sp.symbols("eps", real=True)
eta = sp.diag(-1, 1, 1, 1)

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])

P = sp.Matrix(sp.symbols("P0:4", real=True))
p0 = sp.Matrix(sp.symbols("z0:4", real=True))

def dot(X, Y):
    return (X.T*eta*Y)[0]

def geometry(uu, qq):
    speed = sp.sqrt(-dot(uu, uu))
    T = uu/speed
    w = qq+T*dot(T, qq)
    N = w/sp.sqrt(dot(w, w))
    kap = sp.sqrt(dot(w, w))/speed**2
    return speed, T, N, kap

speed, T, N, kap = geometry(u, q)

assert sp.simplify(speed-v) == 0
assert sp.simplify(kap-k) == 0

print("PASS: adapted geometry", flush=True)

# Canonical equations:
# udot = dHT/dp1 = q
# qdot = dHT/dP
#
# Multipliers are external parameters and are
# held fixed during canonical differentiation.

udot = q

qdot = (
    -v**5*k**2*eta*P/(2*C)
    +LT*T+LN*N
)

# Compute the entire normal-vector derivative
# before restricting the phase-space momenta.

N_eps = geometry(
    u+eps*udot,
    q+eps*qdot
)[2]

Ndot = sp.Matrix([
    sp.diff(N_eps[i], eps).subs(eps, 0)
    for i in range(4)
])

Ndot = sp.Matrix([
    sp.factor(component)
    for component in Ndot
])

print("PASS: all four Ndot components constructed",
      flush=True)

# H_T is independent of y, hence p0dot=0.
# n=p0.N therefore has canonical derivative
# ndot=p0.Ndot.

ndot_full = (p0.T*Ndot)[0]

print("PASS: full canonical ndot assembled",
      flush=True)

r, s, b, c, n = sp.symbols(
    "r s b c n", real=True
)

R = r**2+s**2

p0_time = (
    A*k**2
    +v**4*k**2*R/(4*C)
    -B
)

surface = {
    P[0]: 0,
    P[1]: 0,
    P[2]: r,
    P[3]: s,
    p0[0]: p0_time,
    p0[1]: n,
    p0[2]: b,
    p0[3]: c
}

# Apply constraints after the canonical
# normal-vector derivative is constructed.

actual = sp.cancel(
    sp.together(ndot_full.subs(surface))
)

expected = (
    v*k*p0_time
    -v**3*k*(b*r+c*s)/(2*C)
)

residual = sp.cancel(
    sp.together(actual-expected)
)

print("CANONICAL_NDOT =", actual, flush=True)
print("T74_RESIDUAL =", residual, flush=True)

if residual != 0:
    raise AssertionError(
        "STOP: canonical ndot differs from T74"
    )

print("PASS: exact T74 ndot agreement", flush=True)

# Independent geometric checks:
# N.N=1 implies N.Ndot=0.
# T.N=0 implies
# Tdot.N + T.Ndot = 0.

Tdot = sp.Matrix([
    sp.diff(
        geometry(u+eps*udot, q+eps*qdot)[1][i],
        eps
    ).subs(eps, 0)
    for i in range(4)
])

if sp.simplify(dot(N, Ndot)) != 0:
    raise AssertionError(
        "STOP: N normalization derivative"
    )

if sp.simplify(
    dot(Tdot, N)+dot(T, Ndot)
) != 0:
    raise AssertionError(
        "STOP: tangent-normal orthogonality derivative"
    )

print("PASS: normal normalization preserved")
print("PASS: tangent-normal orthogonality preserved")

# Check the transverse p0 contributions.

if sp.simplify(
    sp.diff(actual, b)+v**3*k*r/(2*C)
) != 0:
    raise AssertionError(
        "STOP: transverse b coefficient"
    )

if sp.simplify(
    sp.diff(actual, c)+v**3*k*s/(2*C)
) != 0:
    raise AssertionError(
        "STOP: transverse c coefficient"
    )

print("PASS: transverse p0 contributions")
print("PASS: constraints imposed after differentiation")
print("PASS: generic A, B, C retained")
print("N TRANSPORT: EXACTLY VERIFIED")
print("INTERNAL TRANSPORT: OPEN")
print("GENERAL E_DOT CANONICAL PROOF: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T83F COMPLETE")
