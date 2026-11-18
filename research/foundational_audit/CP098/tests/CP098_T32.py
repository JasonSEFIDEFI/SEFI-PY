"""CP098-T32: Independent tangential tertiary-bracket audit.

Calculates {psi_T, chi_N} using canonical derivatives
on the general regular external constrained family.

A, B, C are held constant in the internal sector.
Tests an earlier necessary closure identity.
Does not assume that the identity is true.
"""

import sympy as sp

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
a = sp.symbols("a", real=True)
A, B, C = sp.symbols("A B C", real=True, nonzero=True)
r, s, x, y, n, b, c = sp.symbols(
    "r s x y n b c", real=True
)
e = sp.symbols("e", real=True)

def dot(X, Y):
    return (X.T * eta * Y)[0]

def pair(p, X):
    return (p.T * X)[0]

def geom(u, q):
    speed = sp.sqrt(-dot(u, u))
    T = u/speed
    tq = dot(T, q)
    w = q + T*tq
    w2 = dot(w, w)
    kap = sp.sqrt(w2)/speed**2
    N = w/sp.sqrt(w2)
    return speed, T, tq, N, kap, w2

def chi(u, q, p1):
    speed, T, tq, N, kap, w2 = geom(u, q)
    return pair(p1, N) + 2*A*kap/speed

def psi(u, q, P, p0):
    speed, T, tq, N, kap, w2 = geom(u, q)
    return (
        pair(p0, T)
        + B - A*kap**2
        - w2*dot(P, P)/(4*C)
    )

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])
P = sp.Matrix([0, 0, r, s])
p1 = sp.Matrix([0, -2*A*k/v, x, y])

R = r**2 + s**2

p0 = sp.Matrix([
    A*k**2 + v**4*k**2*R/(4*C) - B,
    n, b, c
])

speed, T, tq, N, kap, w2 = geom(u, q)

print("=== CP098-T32 ===", flush=True)

checks = {
    "phi_T": pair(P, T),
    "phi_N": pair(P, N),
    "chi_T": pair(p1, T),
    "chi_N": chi(u, q, p1),
    "psi_T": psi(u, q, P, p0),
}

for name, expr in checks.items():
    result = sp.factor(expr)
    print(name, "=", result, flush=True)
    assert result == 0, name

print("PASS: five constraints", flush=True)

# Canonical pairs:
# (u,p1) and (q,P).
#
# psi does not depend on p1.
# chi does not depend on P.
#
# Therefore:
#
# {psi,chi} =
#     psi_u . chi_p1 - psi_P . chi_q
#
# chi_p1 = N.
#
# All momenta remain fixed during differentiation.

print("Computing psi_u . N...", flush=True)

term_u = sp.factor(
    sp.diff(
        psi(u + e*N, q, P, p0), e
    ).subs(e, 0)
)

print("Computing -psi_P . chi_q...", flush=True)

# psi_P = -w2/(2*C) * eta*P.
# Thus -psi_P . chi_q is the positive term below.
term_q = sp.factor(
    w2/(2*C) *
    sp.diff(
        chi(u, q + e*(eta*P), p1), e
    ).subs(e, 0)
)

direct = sp.factor(term_u + term_q)

print("Computing independently reconstructed F_N...",
      flush=True)

# External F_N from T25-T30, with D_A=0.
FN = (
    -pair(p0, N)
    + 6*A*kap*tq/speed**2
    - speed**3*kap*dot(p1, eta*P)/(2*C)
    + speed**2*kap*tq*dot(P, P)/(2*C)
)

# dot(p1, eta*P) equals p1.T*P because
# eta*eta = identity. This is NOT the desired
# covector-covector Minkowski contraction.
#
# Use the correct covector-covector contraction.
FN = (
    -pair(p0, N)
    + 6*A*kap*tq/speed**2
    - speed**3*kap*dot(p1, P)/(2*C)
    + speed**2*kap*tq*dot(P, P)/(2*C)
)

candidate = sp.factor(
    4*A*kap*tq/speed**3 - FN/speed
)

difference = sp.factor(direct - candidate)

print("psi_u.N =", term_u, flush=True)
print("-psi_P.chi_q =", term_q, flush=True)
print("DIRECT =", direct, flush=True)
print("F_N =", sp.factor(FN), flush=True)
print("CANDIDATE =", candidate, flush=True)
print("DIFFERENCE =", difference, flush=True)

if difference != 0:
    print("STOP: tertiary identity not verified")
    print("Investigate before any closure claim")
    raise SystemExit(1)

print("PASS: external tertiary-bracket identity")
print("INTERNAL SECTOR: NOT YET INCLUDED")
print("GENERAL DIRAC CLOSURE: NOT ESTABLISHED")
print("DYNAMICAL CONSISTENCY: OPEN")
