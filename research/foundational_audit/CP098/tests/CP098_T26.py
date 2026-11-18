"""CP098-T26: External normal-preservation bracket.

Symbolic A, C, B; C nonzero.
Nine-parameter adapted geometry/momentum family.

Computes external canonical contributions independently.
The proposed F_N is introduced only for comparison.
No unrestricted Dirac closure claim.
"""

import sympy as sp

v, k = sp.symbols("v k", positive=True)
A, B, C = sp.symbols("A B C", real=True, nonzero=True)
a, n, r, s, x, y = sp.symbols(
    "a n r s x y", real=True
)
e, z2, z3 = sp.symbols("e z2 z3", real=True)

eta = sp.diag(-1, 1, 1, 1)

def mdot(X, Y):
    return (X.T * eta * Y)[0]

def contract(p, X):
    return (p.T * X)[0]

def geometry(u, q):
    speed = sp.sqrt(-mdot(u, u))
    T = u/speed
    w = q + T*mdot(T, q)
    w2 = mdot(w, w)
    kap = sp.sqrt(w2)/speed**2
    N = w/sp.sqrt(w2)
    return speed, T, N, kap, w2

R = r**2 + s**2
cross = r*x + s*y

u0 = sp.Matrix([v, 0, 0, 0])
q0 = sp.Matrix([a, v**2*k, 0, 0])
P0 = sp.Matrix([0, 0, r, s])
p10 = sp.Matrix([0, -2*A*k/v, x, y])

p00 = A*k**2 + v**4*k**2*R/(4*C) - B
p0 = sp.Matrix([p00, n, 0, 0])

def chi(u, q):
    speed, T, N, kap, w2 = geometry(u, q)
    return contract(p10, N) + 2*A*kap/speed

def H(u, q):
    speed, T, N, kap, w2 = geometry(u, q)
    return (
        contract(p0, u)
        + contract(p10, q)
        + speed*(B + A*kap**2)
        - speed*w2*mdot(P0, P0)/(4*C)
    )

print("=== CP098-T26 ===", flush=True)

speed, T, N, kap, w2 = geometry(u0, q0)

constraints = {
    "phi_T": contract(P0, T),
    "phi_N": contract(P0, N),
    "chi_T": contract(p10, T),
    "chi_N": chi(u0, q0),
    "psi_T": (
        contract(p0, T)
        + B
        - A*kap**2
        - w2*mdot(P0, P0)/(4*C)
    ),
}

for name, expression in constraints.items():
    result = sp.factor(expression)
    print(name, "=", result, flush=True)
    assert result == 0, name

# (u,p1): chi_u.H_p1 - chi_p1.H_u
# H_p1=q and chi_p1=N.

print("Computing u,p1...", flush=True)

du_chi = sp.diff(
    chi(u0 + e*q0, q0), e
).subs(e, 0)

du_H = sp.diff(
    H(u0 + e*N, q0), e
).subs(e, 0)

part_u = sp.factor(du_chi - du_H)
print("(u,p1) =", part_u, flush=True)

# (q,P): chi_q.H_P
# chi_P is identically zero.
# At this constrained point, H_P0=H_P1=0.

print("Computing q,P...", flush=True)

q_shift = sp.Matrix([a, v**2*k, z2, z3])
chi_shift = chi(u0, q_shift)

subs0 = {z2: 0, z3: 0}

dq2 = sp.diff(chi_shift, z2).subs(subs0)
dq3 = sp.diff(chi_shift, z3).subs(subs0)

dH_dP2 = -v**5*k**2*r/(2*C)
dH_dP3 = -v**5*k**2*s/(2*C)

part_q = sp.factor(
    dq2*dH_dP2 + dq3*dH_dP3
)
print("(q,P) =", part_q, flush=True)

direct = sp.factor(part_u + part_q)

# Proposed EXTERNAL terms only.
# Internal contribution 2*k*D_A is excluded.
proposed = (
    -n
    - 6*A*k*a/v**2
    - v**3*k*cross/(2*C)
    - v**2*k*a*R/(2*C)
)

difference = sp.factor(direct - proposed)

print("DIRECT EXTERNAL F_N =", direct, flush=True)
print("PROPOSED EXTERNAL F_N =", sp.factor(proposed), flush=True)
print("DIFFERENCE =", difference, flush=True)

assert difference == 0

print("PASS: symbolic-coupling external family")
print("INTERNAL CONTRIBUTION: separately verified T25")
print("GENERAL NORMAL PRESERVATION: NOT ESTABLISHED")
print("GENERAL DIRAC CLOSURE: NOT ESTABLISHED")
