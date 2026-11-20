"""CP098-T83D: Independent canonical transport of P².

Compute Pdot = -dH_T/dq by directional derivatives.
Then derive {P²,H_T} and restrict to the regular
constraint surface only after differentiation.
"""

import sympy as sp

print("=== CP098-T83D ===", flush=True)

v, k = sp.symbols("v k", positive=True)
a, LN = sp.symbols("a LN", real=True)

A, C = sp.symbols(
    "A C", real=True, nonzero=True
)

eps = sp.symbols("eps", real=True)

eta = sp.diag(-1, 1, 1, 1)

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])

P = sp.Matrix(sp.symbols("P0:4", real=True))
p1 = sp.Matrix(sp.symbols("j0:4", real=True))

basis = sp.eye(4)

def dot(X, Y):
    return (X.T*eta*Y)[0]

def geom(uu, qq):
    speed = sp.sqrt(-dot(uu, uu))
    T = uu/speed
    w = qq+T*dot(T, qq)
    kap = sp.sqrt(dot(w, w))/speed**2
    N = w/sp.sqrt(dot(w, w))
    return speed, T, N, kap

# Only the q-dependent terms in H_T matter
# for Pdot. The remaining Hamiltonian terms
# have zero q derivative.
#
# LT*phi_T is q-independent.
# LN is an external multiplier, held fixed.

def H_q(qq):
    speed, T, N, kap = geom(u, qq)

    return (
        (p1.T*qq)[0]
        +speed*A*kap**2
        -speed**5*kap**2*dot(P, P)/(4*C)
        +LN*(P.T*N)[0]
    )

Pdot = sp.Matrix([
    -sp.diff(
        H_q(q+eps*basis[:,i]), eps
    ).subs(eps, 0)
    for i in range(4)
])

print("PASS: canonical Pdot constructed",
      flush=True)

# P² is the full Lorentzian momentum norm.
# No constraints have yet been imposed.

R_full = dot(P, P)

Rdot_full = 2*dot(P, Pdot)

print("PASS: full momentum-norm derivative built",
      flush=True)

# Apply the normal and tangent momentum
# constraints only after canonical differentiation.

r, s, x, y = sp.symbols(
    "r s x y", real=True
)

surface = {
    P[0]: 0,
    P[1]: 0,
    P[2]: r,
    P[3]: s,
    p1[2]: x,
    p1[3]: y
}

actual = sp.factor(
    Rdot_full.subs(surface)
)

R = r**2+s**2
S = r*x+s*y

expected = (
    -2*S
    -2*LN*R/(v**2*k)
)

residual = sp.cancel(
    sp.together(actual-expected)
)

if residual != 0:
    raise AssertionError(
        f"STOP: Rdot residual = {residual}"
    )

print("PASS: exact regular-surface Rdot",
      flush=True)

# Independently verify that the normal
# multiplier coefficient is correct.

coefficient = sp.diff(actual, LN)

expected_coefficient = -2*R/(v**2*k)

if sp.simplify(
    coefficient-expected_coefficient
) != 0:
    raise AssertionError(
        "STOP: normal multiplier coefficient"
    )

print("PASS: normal multiplier coefficient",
      flush=True)

# Check that the transverse p1 terms survive.

if sp.simplify(
    sp.diff(actual, x)+2*r
) != 0:
    raise AssertionError(
        "STOP: transverse x coefficient"
    )

if sp.simplify(
    sp.diff(actual, y)+2*s
) != 0:
    raise AssertionError(
        "STOP: transverse y coefficient"
    )

print("PASS: transverse momentum dependence",
      flush=True)

print("PASS: constraints imposed after differentiation")
print("PASS: generic A and C retained")
print("R TRANSPORT: EXACTLY VERIFIED")
print("REMAINING SCALAR TRANSPORT: OPEN")
print("GENERAL E_DOT CANONICAL PROOF: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T83D COMPLETE")
