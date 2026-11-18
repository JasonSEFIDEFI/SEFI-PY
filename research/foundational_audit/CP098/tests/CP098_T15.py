"""CP098-T15: Normal primary-secondary bracket sign.

Independent canonical calculation in an adapted timelike
frame, with no constraint substitutions before differentiation.
"""

import sympy as sp

v, A = sp.symbols("v A", positive=True)
q0, q1, q2, q3 = sp.symbols("q0:4", real=True)
P0, P1, P2, P3 = sp.symbols("P0:4", real=True)
r0, r1, r2, r3 = sp.symbols("r0:4", real=True)

q = sp.Matrix([q0, q1, q2, q3])
P = sp.Matrix([P0, P1, P2, P3])
p1 = sp.Matrix([r0, r1, r2, r3])

eta = sp.diag(-1, 1, 1, 1)
T = sp.Matrix([1, 0, 0, 0])

def mdot(x, y):
    return (x.T * eta * y)[0]

# u = v*T, so the transverse acceleration is:
w = q + T * mdot(T, q)

w2 = sp.expand(mdot(w, w))
assert w2 == q1**2 + q2**2 + q3**2

kappa = sp.sqrt(w2) / v**2
N = w / sp.sqrt(w2)

# p1 and P are covectors.
phiN = (P.T * N)[0]
chiN = (p1.T * N)[0] + 2*A*kappa/v

def pb_qP(f, g):
    return sum(
        sp.diff(f, q[i])*sp.diff(g, P[i])
        - sp.diff(f, P[i])*sp.diff(g, q[i])
        for i in range(4)
    )

forward = pb_qP(chiN, phiN)
reverse = pb_qP(phiN, chiN)

expected = 2*A/v**3

print("=== CP098-T15 ===", flush=True)

# Evaluate at a regular point with nonzero curvature.
point = {
    q0: 1,
    q1: 2,
    q2: 0,
    q3: 0,
    v: 4,
    A: 3,
    r0: 5,
    r1: 7,
    r2: 11,
    r3: 13,
    P0: 17,
    P1: 19,
    P2: 23,
    P3: 29,
}

f = sp.simplify(forward.subs(point))
r = sp.simplify(reverse.subs(point))
e = sp.simplify(expected.subs(point))

print("{chi_N,phi_N} =", f, flush=True)
print("{phi_N,chi_N} =", r, flush=True)
print("Expected +2A/v^3 =", e, flush=True)

assert f == e
assert r == -e
assert sp.simplify(forward + reverse) == 0

print("PASS: both bracket orientations", flush=True)

# Check symbolic equality on the regular adapted domain.
print("Checking symbolic identity...", flush=True)

assert sp.simplify(forward - expected) == 0
assert sp.simplify(reverse + expected) == 0

print("PASS: symbolic identities", flush=True)
print("Correction register CP098-C01 requires revision.")
print("No Git changes beyond saving this script.")
