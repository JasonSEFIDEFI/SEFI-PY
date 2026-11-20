"""CP098-T50: Normal secondary Hamiltonian evolution.

Direct external canonical directional differentiation
with unrestricted momenta in the regular adapted frame.

Adds the exact internal bracket for the specified
quadratic internal kinetic model.

No off-surface chain identity or new constraint assumed.
"""

import sympy as sp

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
a = sp.symbols("a", real=True)
A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)
DA = sp.symbols("DA", real=True)
e = sp.symbols("e", real=True)

Pz = sp.symbols("P0:4", real=True)
p1z = sp.symbols("p10:14", real=True)
p0z = sp.symbols("p00:04", real=True)

P = sp.Matrix(Pz)
p1 = sp.Matrix(p1z)
p0 = sp.Matrix(p0z)

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])

def dot(X, Y):
    return (X.T * eta * Y)[0]

def pair(p, X):
    return (p.T * X)[0]

def geom(uu, qq):
    speed = sp.sqrt(-dot(uu, uu))
    T = uu/speed
    tq = dot(T, qq)
    w = qq + T*tq
    w2 = dot(w, w)
    kap = sp.sqrt(w2)/speed**2
    N = w/sp.sqrt(w2)
    return speed, T, N, kap, w2

def chiN(uu, qq):
    speed, T, N, kap, w2 = geom(uu, qq)
    return pair(p1, N) + 2*A*kap/speed

def H(uu, qq):
    speed, T, N, kap, w2 = geom(uu, qq)
    return (
        pair(p0, uu)
        + pair(p1, qq)
        + speed*(B + A*kap**2)
        - speed*w2*dot(P, P)/(4*C)
    )

def directional(expr):
    return sp.diff(expr, e).subs(e, 0)

print("=== CP098-T50 ===", flush=True)

speed, T, N, kap, w2 = geom(u, q)

assert sp.simplify(speed-v) == 0
assert sp.simplify(kap-k) == 0
assert sp.simplify(chiN(u,q)-(p1z[1]+2*A*k/v)) == 0

print("PASS: regular adapted geometry", flush=True)

# External canonical bracket:
#
# {chiN,H} =
# chiN_u . H_p1
# - chiN_p1 . H_u
# + chiN_q . H_P
# - chiN_P . H_q
#
# H_p1=q, chiN_p1=N, chiN_P=0.

print("Differentiating u contribution...", flush=True)
term_u = sp.factor(
    directional(chiN(u+e*q, q))
)

print("Differentiating p1 contribution...", flush=True)
term_p1 = sp.factor(
    -directional(H(u+e*N, q))
)

print("Differentiating q contribution...", flush=True)
H_P = -speed*w2*eta*P/(2*C)

term_q = sp.factor(
    directional(chiN(u, q+e*H_P))
)

external = sp.factor(term_u+term_p1+term_q)

# Internal bracket:
# chiN_internal = 2*A*k/v
# H_internal = v*(B + A*k^2)
#
# {A,B}=DA and {A,A}=0.
internal = 2*k*DA

full = sp.factor(external+internal)

print("u contribution =", term_u, flush=True)
print("p1 contribution =", term_p1, flush=True)
print("q contribution =", term_q, flush=True)
print("External bracket =", external, flush=True)
print("Internal bracket =", internal, flush=True)
print("Full bracket =", full, flush=True)

# Independently verified constrained family.
r, s, x, y, n, b, c = sp.symbols(
    "r s x y n b c", real=True
)

surface = {
    Pz[0]: 0,
    Pz[1]: 0,
    Pz[2]: r,
    Pz[3]: s,

    p1z[0]: 0,
    p1z[1]: -2*A*k/v,
    p1z[2]: x,
    p1z[3]: y,

    p0z[0]:
        A*k**2
        + v**4*k**2*(r*r+s*s)/(4*C)
        - B,
    p0z[1]: n,
    p0z[2]: b,
    p0z[3]: c
}

reduced = sp.factor(full.subs(surface))

# T34 independently established this forcing.
FN = (
    -n
    - 6*A*a*k/v**2
    - v**3*k*(r*x+s*y)/(2*C)
    - v**2*k*a*(r*r+s*s)/(2*C)
    + 2*k*DA
)

residual = sp.factor(reduced-FN)

print("Surface bracket =", reduced, flush=True)
print("T34 forcing =", sp.factor(FN), flush=True)
print("T34 residual =", residual, flush=True)

assert residual == 0

print("PASS: independent T34 forcing cross-check")

# Normal primary bracket from T40.
chiN_phiN = 2*A/v**3

lambda_N = -FN/chiN_phiN

total_reduced = sp.factor(
    reduced+lambda_N*chiN_phiN
)

print("Normal multiplier =", sp.factor(lambda_N),
      flush=True)
print("Total surface evolution =", total_reduced,
      flush=True)

assert total_reduced == 0

print("PASS: normal multiplier consistency")
print("OFF-SURFACE CHI_N BRACKET: CALCULATED")
print("NORMAL CHAIN CLASSIFICATION: PENDING")
print("GAUGE GENERATOR: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
