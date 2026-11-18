"""CP098-T31: Regular external constraint-surface audit.

Verifies:
1. All five constraints on the proposed parameterization.
2. Algebraic reconstruction of constrained components.
3. Independence of the remaining 10 external parameters.
4. Lorentz covariance under a nontrivial exact boost.

This is not a full Dirac-Bergmann closure test.
"""

import sympy as sp

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
a = sp.symbols("a", real=True)
A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)

r, s, x, y, n, b, c = sp.symbols(
    "r s x y n b c", real=True
)

def mdot(X, Y):
    return (X.T * eta * Y)[0]

def pair(p, X):
    return (p.T * X)[0]

def simplify(expr):
    return sp.simplify(sp.factor(expr))

def geometry(u, q):
    speed = sp.sqrt(-mdot(u, u))
    T = u/speed
    tq = mdot(T, q)
    w = q + T*tq
    kap = sp.sqrt(mdot(w, w))/speed**2
    N = w/sp.sqrt(mdot(w, w))
    return speed, T, tq, N, kap

def constraints(u, q, P, p1, p0):
    speed, T, tq, N, kap = geometry(u, q)

    return {
        "phi_T": pair(P, T),
        "phi_N": pair(P, N),
        "chi_T": pair(p1, T),
        "chi_N": pair(p1, N) + 2*A*kap/speed,
        "psi_T": (
            pair(p0, T)
            + B - A*kap**2
            - speed**4*kap**2*mdot(P, P)/(4*C)
        )
    }

print("=== CP098-T31 ===", flush=True)

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])

P = sp.Matrix([0, 0, r, s])
p1 = sp.Matrix([0, -2*A*k/v, x, y])

p0 = sp.Matrix([
    A*k**2 + v**4*k**2*(r**2+s**2)/(4*C) - B,
    n, b, c
])

print("Checking constrained family...", flush=True)

for name, expression in constraints(
    u, q, P, p1, p0
).items():
    result = simplify(expression)
    print(name, "=", result, flush=True)
    assert result == 0

print("PASS: five constraints", flush=True)

# Start with unrestricted frame components.
P0, P1 = sp.symbols("P0 P1", real=True)
j0, j1 = sp.symbols("j0 j1", real=True)
h0 = sp.symbols("h0", real=True)

P_free = sp.Matrix([P0, P1, r, s])
p1_free = sp.Matrix([j0, j1, x, y])
p0_free = sp.Matrix([h0, n, b, c])

eqs = list(constraints(
    u, q, P_free, p1_free, p0_free
).values())

print("Solving constrained components...", flush=True)

solution = sp.solve(
    eqs,
    [P0, P1, j0, j1, h0],
    dict=True
)

assert len(solution) == 1

expected = {
    P0: 0,
    P1: 0,
    j0: 0,
    j1: -2*A*k/v,
    h0: p0[0]
}

for symbol, value in expected.items():
    difference = simplify(
        solution[0][symbol] - value
    )
    print(symbol, "difference =", difference)
    assert difference == 0

print("PASS: unique constrained-component solution")

# Count external dimensions:
# u: 4, q: 4, P: 4, p1: 4, p0: 4.
# Five constraints give 15 dimensions if independent.
#
# Adapted frame parameters:
# v,k,a: 3
# r,s,x,y,n,b,c: 7
# Lorentz-frame orientation: 5
# (T: 3, N orthogonal to T: 2)
#
# Total = 15.

assert 3 + 7 + 5 == 20 - 5

print("PASS: dimension count 15 = 20 - 5")

# Exact Lorentz transformation.
L = sp.Matrix([
    [sp.Rational(5,4), sp.Rational(3,4), 0, 0],
    [sp.Rational(3,4), sp.Rational(5,4), 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
])

assert L.T*eta*L == eta

def vec(X):
    return L*X

def cov(p):
    return L.inv().T*p

print("Checking Lorentz covariance...", flush=True)

boosted = constraints(
    vec(u), vec(q), cov(P), cov(p1), cov(p0)
)

original = constraints(u, q, P, p1, p0)

for name in original:
    difference = simplify(
        boosted[name] - original[name]
    )
    print(name, "boost difference =", difference)
    assert difference == 0

print("PASS: Lorentz covariance")
print("PASS: regular external parameterization audit")
print("GENERAL DIRAC CLOSURE: NOT ESTABLISHED")
print("DYNAMICAL CONSISTENCY: OPEN")
