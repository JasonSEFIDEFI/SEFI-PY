"""CP098-T40: Primary rows of the constraint matrix.

Independent external canonical directional derivatives.
Compute six primary-versus-secondary/tertiary brackets.

Regular adapted family, arbitrary nonzero A,C.
No rank or gauge classification assumed.
"""
import sympy as sp

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
a, n, b, c, r, s, x, y = sp.symbols(
    "a n b c r s x y", real=True
)
A, B, C = sp.symbols("A B C", nonzero=True, real=True)
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
    return speed, T, N, kap, w2

def constraints(u, q, P, p1, p0):
    speed, T, N, kap, w2 = geom(u, q)
    return [
        pair(P, T),
        pair(P, N),
        pair(p1, T),
        pair(p1, N) + 2*A*kap/speed,
        pair(p0, T) + B - A*kap**2
        - w2*dot(P, P)/(4*C)
    ]

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])
P = sp.Matrix([0, 0, r, s])
p1 = sp.Matrix([0, -2*A*k/v, x, y])
R = r**2 + s**2
p0 = sp.Matrix([
    A*k**2 + v**4*k**2*R/(4*C) - B,
    n, b, c
])

speed, T, N, kap, w2 = geom(u, q)
F = constraints(u, q, P, p1, p0)

def derivative(f, du=None, dq=None):
    uu = u if du is None else u + e*du
    qq = q if dq is None else q + e*dq
    expression = constraints(uu, qq, P, p1, p0)[f]
    return sp.diff(expression, e).subs(e, 0)

# Exact momentum gradients at the evaluation point.
# Covector components remain canonical coordinates.
gp1 = [
    sp.zeros(4,1),
    sp.zeros(4,1),
    T,
    N,
    sp.zeros(4,1)
]

gP = [
    T,
    N,
    sp.zeros(4,1),
    sp.zeros(4,1),
    -w2*eta*P/(2*C)
]

def bracket(i, j):
    # {Fi,Fj} =
    # Fi_u.Fj_p1 - Fi_p1.Fj_u
    # + Fi_q.Fj_P - Fi_P.Fj_q
    terms = [
        derivative(i, du=gp1[j]),
        -derivative(j, du=gp1[i]),
        derivative(i, dq=gP[j]),
        -derivative(j, dq=gP[i])
    ]
    return sp.factor(sum(terms))

names = ["phi_T", "phi_N", "chi_T", "chi_N", "psi_T"]

print("=== CP098-T40 ===", flush=True)

for name, expression in zip(names, F):
    residual = sp.factor(expression)
    print(name, "=", residual, flush=True)
    assert residual == 0

print("PASS: all five surface constraints", flush=True)

results = {}

for i in (0, 1):
    for j in (2, 3, 4):
        print(
            "Computing", names[i], names[j], "...",
            flush=True
        )
        value = bracket(i, j)
        results[(i,j)] = value
        print(
            "{" + names[i] + "," + names[j] + "} =",
            value, flush=True
        )

# Verify two brackets already independently
# established in earlier preservation tests.
assert sp.factor(results[(0,2)]) == 0
assert sp.factor(results[(1,2)]) == 0

# From T34: {chi_N,phi_N}=+2A/v^3.
assert sp.factor(
    results[(1,3)] + 2*A/v**3
) == 0

# From T35: {psi_T,phi_N}=-2A*k/v^2.
assert sp.factor(
    results[(1,4)] - 2*A*k/v**2
) == 0

print("PASS: cross-checks against earlier tests")
print("PRIMARY MATRIX ROWS: CALCULATED")
print("REMAINING BRACKETS: OPEN")
print("MATRIX RANK: NOT YET CALCULATED")
print("GAUGE CLASSIFICATION: OPEN")
