"""CP098-T28: Boosted external (u,p1) bracket.

Independent canonical directional derivatives on a
Lorentz-boosted constrained family.

Symbolic A, B, C, n, r, s, x, y.
Internal coordinates held fixed.

Restricted-family verification, not general closure.
"""

import sympy as sp

A, B, C = sp.symbols("A B C", real=True, nonzero=True)
n, r, s, x, y = sp.symbols("n r s x y", real=True)
e = sp.symbols("e", real=True)

eta = sp.diag(-1, 1, 1, 1)

def mdot(X, Y):
    return (X.T * eta * Y)[0]

def contract(p, X):
    return (p.T * X)[0]

def geometry(u, q):
    v = sp.sqrt(-mdot(u, u))
    T = u/v
    Tq = mdot(T, q)
    w = q + T*Tq
    w2 = mdot(w, w)
    k = sp.sqrt(w2)/v**2
    N = w/sp.sqrt(w2)
    return v, T, Tq, N, k, w2

def chi(u, q, p1):
    v, T, Tq, N, k, w2 = geometry(u, q)
    return contract(p1, N) + 2*A*k/v

def H(u, q, p0, p1, P):
    v, T, Tq, N, k, w2 = geometry(u, q)
    return (
        contract(p0, u)
        + contract(p1, q)
        + v*(B + A*k**2)
        - v*w2*mdot(P, P)/(4*C)
    )

# Exact boost with gamma=5/4 and beta=3/5.
L = sp.Matrix([
    [sp.Rational(5,4), sp.Rational(3,4), 0, 0],
    [sp.Rational(3,4), sp.Rational(5,4), 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
])

# Contravariant vectors transform by L.
# Covectors transform by (L inverse).T.
def boost_vector(X):
    return L*X

def boost_covector(p):
    return L.inv().T*p

R = r**2 + s**2
p00 = A + 64*R/C - B

u = boost_vector(sp.Matrix([4, 0, 0, 0]))
q = boost_vector(sp.Matrix([2, 16, 0, 0]))
P = boost_covector(sp.Matrix([0, 0, r, s]))
p1 = boost_covector(sp.Matrix([0, -A/2, x, y]))
p0 = boost_covector(sp.Matrix([p00, n, 0, 0]))

v, T, Tq, N, k, w2 = geometry(u, q)

print("=== CP098-T28 ===", flush=True)

assert sp.simplify(v-4) == 0
assert sp.simplify(k-1) == 0

constraints = {
    "phi_T": contract(P, T),
    "phi_N": contract(P, N),
    "chi_T": contract(p1, T),
    "chi_N": chi(u, q, p1),
    "psi_T": (
        contract(p0, T)
        + B - A*k**2
        - w2*mdot(P, P)/(4*C)
    ),
}

for name, expr in constraints.items():
    value = sp.factor(expr)
    print(name, "=", value, flush=True)
    assert value == 0, name

print("PASS: all five constraints", flush=True)

# H_p1=q, chi_p1=N.
# Both derivatives are evaluated before setting e=0.

print("Computing chi_u . q...", flush=True)

first = sp.factor(
    sp.diff(chi(u + e*q, q, p1), e).subs(e, 0)
)

print("Computing N . H_u...", flush=True)

second = sp.factor(
    sp.diff(H(u + e*N, q, p0, p1, P), e).subs(e, 0)
)

direct = sp.factor(first-second)

# Candidate u,p1 contribution extracted from the
# reconstructed F_N by removing the q,P and internal
# contributions. Here D_A=0 for external differentiation.
candidate = (
    -contract(p0, N)
    + 6*A*k*Tq/v**2
    + v**2*k*Tq*mdot(P, P)/(2*C)
)

difference = sp.factor(direct-candidate)

print("chi_u . q =", first, flush=True)
print("N . H_u =", second, flush=True)
print("Direct =", direct, flush=True)
print("Candidate =", sp.factor(candidate), flush=True)
print("Difference =", difference, flush=True)

assert difference == 0, "External bracket mismatch"

print("PASS: boosted symbolic external (u,p1) family")
print("UNRESTRICTED IDENTITY: NOT ESTABLISHED")
print("GENERAL DIRAC CLOSURE: NOT ESTABLISHED")
