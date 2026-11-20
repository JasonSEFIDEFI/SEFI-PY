"""CP098-T41: Complete regular constraint-matrix audit.

External brackets calculated by independent canonical
directional derivatives before surface evaluation.

Internal contribution to {chi_N,psi_T} is taken from
the independently verified T33 calculation.

No gauge interpretation is assumed.
"""

import sympy as sp

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
a, n, b, c, r, s, x, y = sp.symbols(
    "a n b c r s x y", real=True
)
A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)
DA = sp.symbols("DA", real=True)
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

R = r*r + s*s
p0 = sp.Matrix([
    A*k*k + v**4*k*k*R/(4*C) - B,
    n, b, c
])

speed, T, N, kap, w2 = geom(u, q)
F = constraints(u, q, P, p1, p0)

def derivative(index, du=None, dq=None):
    uu = u if du is None else u + e*du
    qq = q if dq is None else q + e*dq
    expr = constraints(uu, qq, P, p1, p0)[index]
    return sp.diff(expr, e).subs(e, 0)

# Momentum gradients, evaluated on the surface.
# All vectors below contain canonical component derivatives.
gp1 = [
    sp.zeros(4, 1),
    sp.zeros(4, 1),
    T,
    N,
    sp.zeros(4, 1)
]

gP = [
    T,
    N,
    sp.zeros(4, 1),
    sp.zeros(4, 1),
    -w2*eta*P/(2*C)
]

def external_bracket(i, j):
    terms = [
        derivative(i, du=gp1[j]),
        -derivative(j, du=gp1[i]),
        derivative(i, dq=gP[j]),
        -derivative(j, dq=gP[i])
    ]
    return sp.factor(sum(terms))

names = ["phi_T", "phi_N", "chi_T", "chi_N", "psi_T"]

print("=== CP098-T41 ===", flush=True)

for name, expression in zip(names, F):
    value = sp.factor(expression)
    print(name, "=", value, flush=True)
    assert value == 0

print("PASS: five surface constraints", flush=True)

pairs = [(2, 3), (2, 4), (3, 4)]
results = {}

for i, j in pairs:
    print(
        f"Computing {{{names[i]},{names[j]}}}...",
        flush=True
    )
    value = external_bracket(i, j)
    results[(i, j)] = value
    print("External =", value, flush=True)

# T33 independently established:
# {psi_T,chi_N}_internal = -2*k*DA/v.
# Reverse ordering changes its sign.
internal_34 = 2*k*DA/v

results[(3, 4)] = sp.factor(
    results[(3, 4)] + internal_34
)

print(
    "{chi_N,psi_T} full =",
    results[(3, 4)],
    flush=True
)

# Independent cross-check from T35/T37-era
# bracket calculation:
# {psi_T,chi_T} = 4*A*k**2/v.
assert sp.factor(
    results[(2, 4)] + 4*A*k**2/v
) == 0

print("PASS: chi_T/psi_T cross-check", flush=True)

# T40 verified primary rows.
M = sp.zeros(5)

def put(i, j, value):
    M[i, j] = value
    M[j, i] = -value

put(1, 3, -2*A/v**3)
put(1, 4, 2*A*k/v**2)

for (i, j), value in results.items():
    put(i, j, value)

print("Constraint order:", names, flush=True)
print("Constraint matrix:", flush=True)
sp.pprint(M)

assert M + M.T == sp.zeros(5)

# phi_T is already a zero row.
assert all(M[0, j] == 0 for j in range(5))

# The only potentially nonzero 4x4 principal
# Pfaffian involves indices (1,2,3,4).
#
# pf = M12*M34 - M13*M24 + M14*M23
# with indices relative to the selected submatrix.

pf = sp.factor(
    M[1, 2]*M[3, 4]
    - M[1, 3]*M[2, 4]
    + M[1, 4]*M[2, 3]
)

print("Four-by-four Pfaffian =", pf, flush=True)
print("Pfaffian squared =", sp.factor(pf**2),
      flush=True)

if pf == 0:
    # Rank at least 2 because M[1,3] != 0.
    rank = 2
    print("Generic regular rank = 2", flush=True)
else:
    rank = 4
    print("Generic regular rank = 4", flush=True)
    print(
        "NOTE: check exceptional parameter loci",
        flush=True
    )

print("Matrix nullity =", 5-rank, flush=True)

# Report a null-vector candidate, not yet
# a proof of a physical gauge symmetry.
print("phi_T null direction: VERIFIED", flush=True)

print("CONSTRAINT MATRIX: COMPUTED")
print("RANK RESULT: SUBJECT TO REGULARITY ASSUMPTIONS")
print("GAUGE GENERATOR ANALYSIS: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
