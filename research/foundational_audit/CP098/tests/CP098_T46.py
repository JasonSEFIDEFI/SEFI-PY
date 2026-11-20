"""CP098-T46: Off-surface primary Hamiltonian brackets.

Canonical directional differentiation with unrestricted
external momenta in the regular adapted geometric frame.

Computes {phi_T,Hc} and {phi_N,Hc}.
No off-surface chain identity is presumed.
"""

import sympy as sp

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
a = sp.symbols("a", real=True)
A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)
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

def geometry(uu, qq):
    speed = sp.sqrt(-dot(uu, uu))
    T = uu/speed
    tq = dot(T, qq)
    w = qq + T*tq
    w2 = dot(w, w)
    N = w/sp.sqrt(w2)
    kap = sp.sqrt(w2)/speed**2
    return speed, T, N, kap, w2

def primary(index, uu, qq):
    speed, T, N, kap, w2 = geometry(uu, qq)
    return pair(P, T if index == 0 else N)

def H(uu, qq):
    speed, T, N, kap, w2 = geometry(uu, qq)
    return (
        pair(p0, uu)
        + pair(p1, qq)
        + speed*(B + A*kap**2)
        - speed*w2*dot(P, P)/(4*C)
    )

def directional(expr):
    return sp.diff(expr, e).subs(e, 0)

speed, T, N, kap, w2 = geometry(u, q)

def primary_H_bracket(index):
    normal = T if index == 0 else N

    # Canonical pairs (u,p1) and (q,P).
    # H_p1 = q.
    term_u = directional(
        primary(index, u + e*q, q)
    )

    # H_P = -speed*w2*eta*P/(2*C).
    H_P = -speed*w2*eta*P/(2*C)
    term_q = directional(
        primary(index, u, q + e*H_P)
    )

    # Primary momentum gradient = T or N.
    term_P = -directional(
        H(u, q + e*normal)
    )

    return (
        sp.factor(term_u),
        sp.factor(term_q),
        sp.factor(term_P)
    )

print("=== CP098-T46 ===", flush=True)

print("Momentum variables unrestricted:", flush=True)
print("P =", P.T)
print("p1 =", p1.T)
print("p0 =", p0.T)

assert sp.simplify(primary(0, u, q) - Pz[0]) == 0
assert sp.simplify(primary(1, u, q) - Pz[1]) == 0

print("PASS: primary definitions", flush=True)

results = []

for index, name in enumerate(["phi_T", "phi_N"]):
    print("Computing {" + name + ",Hc}...",
          flush=True)

    terms = primary_H_bracket(index)
    value = sp.factor(sum(terms))

    print("u contribution =", terms[0], flush=True)
    print("q contribution =", terms[1], flush=True)
    print("P contribution =", terms[2], flush=True)
    print("Full bracket =", value, flush=True)

    results.append(value)

# Constrained-family substitution for the
# independent T38R consistency check.
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

for name, expression in zip(
    ["phi_T", "phi_N"], results
):
    residual = sp.factor(expression.subs(surface))
    print(
        name + " surface residual =",
        residual, flush=True
    )
    assert residual == 0, name

print("PASS: independent T38R surface cross-check")

# Check whether the two off-surface brackets
# are identically zero. This is diagnostic,
# not an assumed requirement.
for name, expression in zip(
    ["phi_T", "phi_N"], results
):
    print(
        name + " identically zero off surface:",
        expression == 0,
        flush=True
    )

print("OFF-SURFACE PRIMARY BRACKETS: CALCULATED")
print("SECONDARY CHAIN IDENTITIES: NOT YET CLASSIFIED")
print("GAUGE GENERATOR: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
