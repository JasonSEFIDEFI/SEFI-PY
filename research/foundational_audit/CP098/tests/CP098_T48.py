"""CP098-T48: Off-surface tangential secondary evolution.

Direct canonical directional differentiation of
{chi_T,Hc}, retaining unrestricted external momenta.

No tertiary identity or gauge structure is assumed.
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

def geom(uu, qq):
    speed = sp.sqrt(-dot(uu, uu))
    T = uu/speed
    tq = dot(T, qq)
    w = qq + T*tq
    w2 = dot(w, w)
    kap = sp.sqrt(w2)/speed**2
    return speed, T, kap, w2

def chiT(uu, qq):
    speed, T, kap, w2 = geom(uu, qq)
    return pair(p1, T)

def H(uu, qq):
    speed, T, kap, w2 = geom(uu, qq)
    return (
        pair(p0, uu)
        + pair(p1, qq)
        + speed*(B + A*kap**2)
        - speed*w2*dot(P, P)/(4*C)
    )

def direction(expr):
    return sp.diff(expr, e).subs(e, 0)

print("=== CP098-T48 ===", flush=True)

speed, T, kap, w2 = geom(u, q)

assert sp.simplify(chiT(u, q)-p1z[0]) == 0
print("PASS: adapted chi_T definition", flush=True)

# Canonical bracket:
#
# {chi_T,Hc}
# = chiT_u.H_p1 - chiT_p1.H_u
# + chiT_q.H_P - chiT_P.H_q
#
# chiT_q = chiT_P = 0.
# H_p1 = q and chiT_p1 = T.

term_u = sp.factor(
    direction(chiT(u+e*q, q))
)

term_p1 = sp.factor(
    -direction(H(u+e*T, q))
)

full = sp.factor(term_u+term_p1)

print("u contribution =", term_u, flush=True)
print("p1 contribution =", term_p1, flush=True)
print("{chi_T,Hc} =", full, flush=True)

# Primary + secondary constraints:
#
# P0=P1=0
# p10=0
# p11=-2*A*k/v
#
# No tertiary condition on p00 is imposed.

primary_secondary = {
    Pz[0]: 0,
    Pz[1]: 0,
    p1z[0]: 0,
    p1z[1]: -2*A*k/v
}

reduced = sp.factor(
    full.subs(primary_secondary)
)

print(
    "Primary-secondary surface residual =",
    reduced,
    flush=True
)

# The T37 result predicts preservation on the
# complete five-constraint surface.
r, s, x, y, n, b, c = sp.symbols(
    "r s x y n b c", real=True
)

complete_surface = dict(primary_secondary)
complete_surface.update({
    Pz[2]: r,
    Pz[3]: s,
    p1z[2]: x,
    p1z[3]: y,
    p0z[0]:
        A*k**2
        + v**4*k**2*(r*r+s*s)/(4*C)
        - B,
    p0z[1]: n,
    p0z[2]: b,
    p0z[3]: c
})

final_residual = sp.factor(
    full.subs(complete_surface)
)

print(
    "Complete-surface residual =",
    final_residual,
    flush=True
)

assert final_residual == 0

print("PASS: independent T37 surface cross-check")

# Test whether the tertiary constraint is needed
# for this particular preservation equation.
if reduced == 0:
    print(
        "RESULT: chi_T preserved already on "
        "primary-secondary surface",
        flush=True
    )
else:
    print(
        "RESULT: chi_T requires further "
        "constraint analysis",
        flush=True
    )

print("OFF-SURFACE CHI_T BRACKET: CALCULATED")
print("CHI_N OFF-SURFACE EVOLUTION: OPEN")
print("GAUGE GENERATOR: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
