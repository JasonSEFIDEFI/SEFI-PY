"""CP098-T35: Direct external tangential-preservation brackets.

Canonical pairs: (u,p1), (q,P).
External sector only; internal variables held fixed.
All derivatives are taken before constraint substitution.
No assumption that tertiary consistency closes.
"""

import sympy as sp

eta = sp.diag(-1, 1, 1, 1)
v, k = sp.symbols("v k", positive=True)
a, n, b, c, r, s, x, y = sp.symbols(
    "a n b c r s x y", real=True
)
A, B, C = sp.symbols("A B C", real=True, nonzero=True)
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

def psi(u, q, P, p0):
    speed, T, tq, N, kap, w2 = geom(u, q)
    return pair(p0, T) + B - A*kap**2 - w2*dot(P, P)/(4*C)

def H(u, q, P, p1, p0):
    speed, T, tq, N, kap, w2 = geom(u, q)
    return (
        pair(p0, u) + pair(p1, q)
        + speed*(B + A*kap**2)
        - speed*w2*dot(P, P)/(4*C)
    )

def phi_T(u, q, P):
    return pair(P, geom(u, q)[1])

def phi_N(u, q, P):
    return pair(P, geom(u, q)[3])

def directional(f):
    return sp.factor(sp.diff(f, e).subs(e, 0))

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])
P = sp.Matrix([0, 0, r, s])
p1 = sp.Matrix([0, -2*A*k/v, x, y])
R = r**2+s**2
p0 = sp.Matrix([
    A*k**2 + v**4*k**2*R/(4*C)-B,
    n, b, c
])

speed, T, tq, N, kap, w2 = geom(u, q)

print("=== CP098-T35 ===", flush=True)

checks = {
    "phi_T": phi_T(u, q, P),
    "phi_N": phi_N(u, q, P),
    "chi_T": pair(p1, T),
    "chi_N": pair(p1, N)+2*A*kap/speed,
    "psi_T": psi(u, q, P, p0),
}

for name, expression in checks.items():
    value = sp.factor(expression)
    print(name, "=", value, flush=True)
    assert value == 0, name

print("PASS: five constraints", flush=True)

# {psi,H}_external:
#
# psi_p1=0, H_p1=q.
# psi_P=-w2*eta*P/(2C).
# H_P=-speed*w2*eta*P/(2C).
#
# Thus three independently differentiated terms:
# psi_u . q
# psi_q . H_P
# -psi_P . H_q

print("Calculating {psi,Hc}_external...", flush=True)

term_u = directional(
    psi(u+e*q, q, P, p0)
)

term_qP = sp.factor(
    -speed*w2/(2*C) *
    directional(psi(u, q+e*(eta*P), P, p0))
)

term_Pq = sp.factor(
    w2/(2*C) *
    directional(H(u, q+e*(eta*P), P, p1, p0))
)

bracket_H = sp.factor(term_u+term_qP+term_Pq)

print("psi_u.H_p1 =", term_u, flush=True)
print("psi_q.H_P =", term_qP, flush=True)
print("-psi_P.H_q =", term_Pq, flush=True)
print("{psi,Hc}_external =", bracket_H, flush=True)

# {psi,phi_T} and {psi,phi_N}:
# psi has no p1 dependence, and phi has no p1
# dependence. Only (q,P) contributes.
#
# phi_P = T or N.
# -psi_P.phi_q = +w2/(2C) directional phi
# along eta*P.

print("Calculating multiplier brackets...", flush=True)

bracket_phiT = sp.factor(
    directional(psi(u, q+e*T, P, p0))
    + w2/(2*C)*directional(
        phi_T(u, q+e*(eta*P), P)
    )
)

bracket_phiN = sp.factor(
    directional(psi(u, q+e*N, P, p0))
    + w2/(2*C)*directional(
        phi_N(u, q+e*(eta*P), P)
    )
)

print("{psi,phi_T} =", bracket_phiT, flush=True)
print("{psi,phi_N} =", bracket_phiN, flush=True)

# These are testable geometric predictions,
# not substitutions used to calculate brackets.
assert sp.factor(bracket_phiT) == 0
assert sp.factor(bracket_phiN + 2*A*k/v**2) == 0

print("PASS: independent multiplier brackets")

# External normal-preservation bracket from T29/T30.
FN_ext = (
    -n - 6*A*a*k/v**2
    - v**3*k*(r*x+s*y)/(2*C)
    - v**2*k*a*R/(2*C)
)

# This lambda_N is ONLY the external-sector value.
# Full lambda_N requires the internal contribution.
lambda_N_ext = -v**3*FN_ext/(2*A)

residual_ext = sp.factor(
    bracket_H + lambda_N_ext*bracket_phiN
)

print("lambda_N (external-only) =", sp.factor(lambda_N_ext),
      flush=True)
print("External consistency residual =", residual_ext,
      flush=True)

print("EXTERNAL BRACKETS: CALCULATED")
print("INTERNAL CONTRIBUTIONS: NOT INCLUDED")
print("TERTIARY CONSISTENCY: NOT YET RESOLVED")
print("DYNAMICAL CONSISTENCY: OPEN")
