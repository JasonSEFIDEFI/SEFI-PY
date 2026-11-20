"""CP098-T43: Functional independence and basis regularity.

Check the constraint Jacobian with respect to five
independent canonical momentum components.

Also check that the proposed 3+2 constraint basis
is invertible on the regular sector.

No gauge or physical-degree count is assumed.
"""

import sympy as sp

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)
a, n, b, c, r, s, x, y = sp.symbols(
    "a n b c r s x y", real=True
)

P0, P1, p10, p11, p00 = sp.symbols(
    "P0 P1 p10 p11 p00", real=True
)

def dot(X, Y):
    return (X.T * eta * Y)[0]

def pair(p, X):
    return (p.T * X)[0]

def constraints(u, q, P, p1, p0):
    speed = sp.sqrt(-dot(u, u))
    T = u/speed
    tq = dot(T, q)
    w = q + T*tq
    w2 = dot(w, w)
    kap = sp.sqrt(w2)/speed**2
    N = w/sp.sqrt(w2)

    return sp.Matrix([
        pair(P, T),
        pair(P, N),
        pair(p1, T),
        pair(p1, N) + 2*A*kap/speed,
        pair(p0, T) + B - A*kap**2
        - w2*dot(P, P)/(4*C)
    ])

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])

P = sp.Matrix([P0, P1, r, s])
p1 = sp.Matrix([p10, p11, x, y])
p0 = sp.Matrix([p00, n, b, c])

F = constraints(u, q, P, p1, p0)

momenta = [P0, P1, p10, p11, p00]

J = F.jacobian(momenta)

surface = {
    P0: 0,
    P1: 0,
    p10: 0,
    p11: -2*A*k/v,
    p00: A*k**2
        + v**4*k**2*(r*r+s*s)/(4*C) - B
}

Js = sp.simplify(J.subs(surface))

print("=== CP098-T43 ===", flush=True)

for name, f in zip(
    ["phi_T", "phi_N", "chi_T", "chi_N", "psi_T"],
    F
):
    residual = sp.simplify(f.subs(surface))
    print(name, "=", residual, flush=True)
    assert residual == 0

print("PASS: constrained family", flush=True)

print("Momentum Jacobian:", flush=True)
sp.pprint(Js)

detJ = sp.factor(Js.det())
print("Jacobian determinant =", detJ, flush=True)

assert detJ != 0
print("PASS: five functionally independent constraints",
      flush=True)

# T42 basis:
#
# G1 = phi_T
# G2 = chi_T + 2*t*phi_N
# G3 = psi_T + t*chi_N - (E/d)*phi_N
#
# S1 = phi_N
# S2 = chi_N

d, t, E = sp.symbols(
    "d t E", real=True, nonzero=True
)

basis = sp.Matrix([
    [1, 0,    0, 0, 0],
    [0, 2*t,  1, 0, 0],
    [0, -E/d, 0, t, 1],
    [0, 1,    0, 0, 0],
    [0, 0,    0, 1, 0]
])

det_basis = sp.factor(basis.det())

print("Basis transformation determinant =",
      det_basis, flush=True)

assert det_basis != 0
print("PASS: invertible 3+2 constraint basis",
      flush=True)

# Verify the bracket matrix in the new basis.
M = sp.Matrix([
    [0, 0, 0, 0, 0],
    [0, 0, 0, -d, t*d],
    [0, 0, 0, 2*t*d, -2*t*t*d],
    [0, d, -2*t*d, 0, E],
    [0, -t*d, 2*t*t*d, -E, 0]
])

newM = sp.simplify(basis*M*basis.T)

expected = sp.zeros(5)
expected[3,4] = -d
expected[4,3] = d

print("Transformed bracket matrix:", flush=True)
sp.pprint(newM)

assert newM == expected

print("PASS: first/second-class block decomposition")
print("FUNCTIONAL INDEPENDENCE: VERIFIED LOCALLY")
print("REGULAR CONSTRAINT BASIS: VERIFIED")
print("GAUGE GENERATOR: OPEN")
print("PHYSICAL DEGREES OF FREEDOM: NOT YET COUNTED")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
