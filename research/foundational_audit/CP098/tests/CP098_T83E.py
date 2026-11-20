"""CP098-T83E: Independent canonical transport of S=P.p1.

All four components of Pdot and p1dot are included.
Canonical differentiation precedes surface substitution.
No T74 code or transport rates are imported.
"""

import sympy as sp

print("=== CP098-T83E ===", flush=True)

v, k = sp.symbols("v k", positive=True)
a, LT, LN = sp.symbols("a LT LN", real=True)

A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)

eps = sp.symbols("eps", real=True)
eta = sp.diag(-1, 1, 1, 1)
basis = sp.eye(4)

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])

P = sp.Matrix(sp.symbols("P0:4", real=True))
p1 = sp.Matrix(sp.symbols("j0:4", real=True))
p0 = sp.Matrix(sp.symbols("z0:4", real=True))

def dot(X, Y):
    return (X.T*eta*Y)[0]

def pair(momentum, X):
    return (momentum.T*X)[0]

def geom(uu, qq):
    speed = sp.sqrt(-dot(uu, uu))
    T = uu/speed
    w = qq+T*dot(T, qq)
    kap = sp.sqrt(dot(w, w))/speed**2
    N = w/sp.sqrt(dot(w, w))
    return speed, T, N, kap

def HT(uu, qq):
    speed, T, N, kap = geom(uu, qq)

    return (
        pair(p0, uu)
        +pair(p1, qq)
        +speed*(B+A*kap**2)
        -speed**5*kap**2*dot(P, P)/(4*C)
        +LT*pair(P, T)
        +LN*pair(P, N)
    )

# Directional canonical derivatives at the
# adapted point, with all momenta unrestricted.

Pdot = []
p1dot = []

for i in range(4):
    component = -sp.diff(
        HT(u, q+eps*basis[:, i]), eps
    ).subs(eps, 0)

    Pdot.append(sp.factor(component))

print("PASS: all four Pdot components built",
      flush=True)

for i in range(4):
    component = -sp.diff(
        HT(u+eps*basis[:, i], q), eps
    ).subs(eps, 0)

    p1dot.append(sp.factor(component))

print("PASS: all four p1dot components built",
      flush=True)

# S=P.p1 uses the ordinary canonical momentum
# contraction, not the Minkowski metric.
#
# Therefore:
# Sdot = Pdot.p1 + P.p1dot

Sdot_full = sum(
    Pdot[i]*p1[i]+P[i]*p1dot[i]
    for i in range(4)
)

print("PASS: full canonical Sdot assembled",
      flush=True)

r, s, x, y, b, c, n = sp.symbols(
    "r s x y b c n", real=True
)

surface = {
    P[0]: 0,
    P[1]: 0,
    P[2]: r,
    P[3]: s,
    p1[0]: 0,
    p1[1]: -2*A*k/v,
    p1[2]: x,
    p1[3]: y,
    p0[1]: n,
    p0[2]: b,
    p0[3]: c
}

# Surface restriction only after full
# canonical differentiation.
actual = sp.cancel(
    sp.together(Sdot_full.subs(surface))
)

R = r**2+s**2
S = r*x+s*y

# Independently supplied T74 target:
expected = (
    -x**2-y**2-r*b-s*c
    -LN*S/(v**2*k)
    -LT*R/v
    +a*LN*R/(k*v**3)
    -A*k**2*v**2*R/C
)

residual = sp.cancel(
    sp.together(actual-expected)
)

print("CANONICAL_SDOT =", actual, flush=True)
print("T74_RESIDUAL =", residual, flush=True)

if residual != 0:
    raise AssertionError(
        "STOP: canonical Sdot differs from T74"
    )

print("PASS: exact T74 Sdot agreement",
      flush=True)

# Check the contribution from the longitudinal
# component Pdot[1]*p1[1] independently.

longitudinal = sp.factor(
    (Pdot[1]*p1[1]).subs(surface)
)

expected_longitudinal = (
    -A*k**2*v**2*R/C
)

if sp.cancel(
    sp.together(
        longitudinal-expected_longitudinal
    )
) != 0:
    raise AssertionError(
        "STOP: longitudinal contribution mismatch"
    )

print("PASS: longitudinal momentum contribution",
      flush=True)

print("PASS: constraints applied after differentiation")
print("PASS: generic A, B, C retained")
print("S TRANSPORT: EXACTLY VERIFIED")
print("REMAINING SCALAR TRANSPORT: OPEN")
print("GENERAL E_DOT CANONICAL PROOF: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T83E COMPLETE")
