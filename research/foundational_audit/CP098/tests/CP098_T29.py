"""CP098-T29: Symbolic regular constrained-family audit.

Tests the external (u,p1) contribution for symbolic
v, kappa, acceleration tangential component, couplings,
and transverse momenta.

Does not establish unrestricted off-surface identity.
"""
import sympy as sp

eta = sp.diag(-1, 1, 1, 1)
v, k, a = sp.symbols("v k a", positive=True)
A, B, C = sp.symbols("A B C", nonzero=True, real=True)
n, r, s, x, y = sp.symbols("n r s x y", real=True)
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

def H(u, q, p0, p1, P):
    speed, T, tq, N, kap, w2 = geom(u, q)
    return (
        pair(p0, u) + pair(p1, q)
        + speed*(B + A*kap**2)
        - speed*w2*dot(P, P)/(4*C)
    )

# General regular adapted family, with a free tangential
# acceleration component. A constant Lorentz boost does
# not change the scalar bracket.
u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])
P = sp.Matrix([0, 0, r, s])
p1 = sp.Matrix([0, -2*A*k/v, x, y])

R = r**2 + s**2
p0 = sp.Matrix([
    A*k**2 + v**4*k**2*R/(4*C) - B,
    n, 0, 0
])

speed, T, tq, N, kap, w2 = geom(u, q)

print("=== CP098-T29 ===", flush=True)

checks = {
    "phi_T": pair(P, T),
    "phi_N": pair(P, N),
    "chi_T": pair(p1, T),
    "chi_N": chi(u, q, p1),
    "psi_T": (
        pair(p0, T) + B - A*kap**2
        - w2*dot(P, P)/(4*C)
    )
}

for name, expression in checks.items():
    result = sp.simplify(expression)
    print(name, "=", result, flush=True)
    assert result == 0, name

print("PASS: five constraints", flush=True)

print("Computing chi_u . q...", flush=True)
first = sp.simplify(
    sp.diff(chi(u+e*q, q, p1), e).subs(e, 0)
)

print("Computing N . H_u...", flush=True)
second = sp.simplify(
    sp.diff(H(u+e*N, q, p0, p1, P), e).subs(e, 0)
)

direct = sp.factor(first-second)

candidate = sp.factor(
    -pair(p0, N)
    + 6*A*kap*tq/speed**2
    + speed**2*kap*tq*dot(P, P)/(2*C)
)

difference = sp.factor(direct-candidate)

print("chi_u . q =", first, flush=True)
print("N . H_u =", second, flush=True)
print("Direct =", direct, flush=True)
print("Candidate =", candidate, flush=True)
print("Difference =", difference, flush=True)

assert difference == 0, "T29 mismatch"

print("PASS: symbolic regular constrained family")
print("UNRESTRICTED IDENTITY: NOT ESTABLISHED")
print("GENERAL DIRAC CLOSURE: NOT ESTABLISHED")
