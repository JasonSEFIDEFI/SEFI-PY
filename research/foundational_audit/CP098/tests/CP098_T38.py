"""CP098-T38: Primary preservation, rebuilt.

Independent canonical variables throughout differentiation.
Substitute the regular constrained family only afterward.

Canonical pairs:
    (u^mu, p1_mu)
    (q^mu, P_mu)

Internal contributions vanish for primary constraints
because phi_T and phi_N have no internal dependence.

No closure assumed.
"""

import sympy as sp

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
a, n, b, c, r, s, x, y = sp.symbols(
    "a n b c r s x y", real=True
)

A, B, C = sp.symbols("A B C", real=True, nonzero=True)
DA = sp.symbols("DA", real=True)
LT, LN = sp.symbols("lambda_T lambda_N", real=True)

u = sp.Matrix(sp.symbols("u0:4", real=True))
q = sp.Matrix(sp.symbols("q0:4", real=True))
P = sp.Matrix(sp.symbols("P0:4", real=True))
p1 = sp.Matrix(sp.symbols("p10:14", real=True))
p0 = sp.Matrix(sp.symbols("p00:04", real=True))

def dot(X, Y):
    return (X.T * eta * Y)[0]

def pair(p, X):
    return (p.T * X)[0]

speed = sp.sqrt(-dot(u, u))
T = u / speed
tq = dot(T, q)
w = q + T*tq
w2 = dot(w, w)
kap = sp.sqrt(w2)/speed**2
N = w/sp.sqrt(w2)

phi_T = pair(P, T)
phi_N = pair(P, N)

Hc = (
    pair(p0, u)
    + pair(p1, q)
    + speed*(B + A*kap**2)
    - speed*w2*dot(P, P)/(4*C)
)

def gradient(f, variables):
    return sp.Matrix([
        sp.diff(f, z) for z in variables
    ])

def PB(f, g):
    # Compute bracket before any surface substitution.
    return (
        pair(gradient(f, u), gradient(g, p1))
        - pair(gradient(f, p1), gradient(g, u))
        + pair(gradient(f, q), gradient(g, P))
        - pair(gradient(f, P), gradient(g, q))
    )

# IMPORTANT:
# Canonical brackets contract component derivatives
# without inserting eta. The pair() helper above
# is a plain component contraction.

R = r**2+s**2

u_surface = [v, 0, 0, 0]
q_surface = [a, v**2*k, 0, 0]
P_surface = [0, 0, r, s]
p1_surface = [0, -2*A*k/v, x, y]
p0_surface = [
    A*k**2 + v**4*k**2*R/(4*C)-B,
    n, b, c
]

values = {}
for symbols, components in [
    (u, u_surface),
    (q, q_surface),
    (P, P_surface),
    (p1, p1_surface),
    (p0, p0_surface),
]:
    values.update(zip(symbols, components))

def surface(expr):
    return sp.factor(expr.subs(values, simultaneous=True))

print("=== CP098-T38 REBUILT ===", flush=True)

print("Checking primary constraints...", flush=True)

assert surface(phi_T) == 0
assert surface(phi_N) == 0

print("PASS: primary constraints", flush=True)

# Compute each bracket independently.
print("Differentiating {phi_T,Hc}...", flush=True)
FT = surface(PB(phi_T, Hc))

print("Differentiating {phi_N,Hc}...", flush=True)
FN = surface(PB(phi_N, Hc))

print("Differentiating primary bracket matrix...", flush=True)
BTT = surface(PB(phi_T, phi_T))
BTN = surface(PB(phi_T, phi_N))
BNT = surface(PB(phi_N, phi_T))
BNN = surface(PB(phi_N, phi_N))

print("{phi_T,Hc} =", FT, flush=True)
print("{phi_N,Hc} =", FN, flush=True)
print("{phi_T,phi_T} =", BTT, flush=True)
print("{phi_T,phi_N} =", BTN, flush=True)
print("{phi_N,phi_T} =", BNT, flush=True)
print("{phi_N,phi_N} =", BNN, flush=True)

assert BTT == 0
assert BNN == 0
assert sp.factor(BTN+BNT) == 0

print("PASS: primary bracket antisymmetry", flush=True)

eqT = sp.factor(FT + LT*BTT + LN*BTN)
eqN = sp.factor(FN + LT*BNT + LN*BNN)

print("dot(phi_T) =", eqT, flush=True)
print("dot(phi_N) =", eqN, flush=True)

# T34 normal multiplier, including internal sector.
FN_chi = (
    -n
    -6*A*a*k/v**2
    -v**3*k*(r*x+s*y)/(2*C)
    -v**2*k*a*R/(2*C)
    +2*k*DA
)

LN_normal = -v**3*FN_chi/(2*A)

after_T = sp.factor(eqT.subs(LN, LN_normal))
after_N = sp.factor(eqN.subs(LN, LN_normal))

print("Normal multiplier =", sp.factor(LN_normal),
      flush=True)
print("dot(phi_T) after multiplier =", after_T,
      flush=True)
print("dot(phi_N) after multiplier =", after_N,
      flush=True)

print("PRIMARY PRESERVATION AUDIT COMPLETED")
print("NO CLOSURE ASSERTION")
print("DYNAMICAL CONSISTENCY: OPEN")
