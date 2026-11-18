"""CP098-T24: Independent normal-preservation audit.

Compute {chi_N,H_c} directly on the nine-parameter
adapted constrained family.

The proposed F_N expression is used ONLY after
the canonical calculation, as a comparison.

Toy internal sector: A=I, C=1, B=ps**2/2.
General Dirac closure is not established.
"""

import sympy as sp

v, k = sp.symbols("v k", positive=True)
a, t, n, r, s, x, y = sp.symbols(
    "a t n r s x y", real=True
)
e, z2, z3 = sp.symbols("e z2 z3", real=True)
I, ps = sp.symbols("I ps", real=True)

eta = sp.diag(-1, 1, 1, 1)

def mdot(X, Y):
    return (X.T * eta * Y)[0]

def contract(p, X):
    return (p.T * X)[0]

def geometry(u, q):
    speed = sp.sqrt(-mdot(u, u))
    T = u / speed
    w = q + T*mdot(T, q)
    w2 = mdot(w, w)
    kap = sp.sqrt(w2)/speed**2
    N = w/sp.sqrt(w2)
    return speed, T, N, kap, w2

u0 = sp.Matrix([v, 0, 0, 0])
q0 = sp.Matrix([a, v**2*k, 0, 0])
P0 = sp.Matrix([0, 0, r, s])
p10 = sp.Matrix([0, -2*k/v, x, y])

R = r**2 + s**2
cross = r*x + s*y

p00 = k**2 + v**4*k**2*R/4 - t**2/2
p0 = sp.Matrix([p00, n, 0, 0])

def chi(u, q):
    speed, T, N, kap, w2 = geometry(u, q)
    return contract(p10, N) + 2*I*kap/speed

def H(u, q):
    speed, T, N, kap, w2 = geometry(u, q)
    return (
        contract(p0, u)
        + contract(p10, q)
        + speed*(ps**2/2 + I*kap**2)
        - speed*w2*mdot(P0, P0)/4
    )

def at_internal(expr):
    return expr.subs({I: 1, ps: t})

print("=== CP098-T24 ===", flush=True)

# Canonical (u,p1) bracket:
# chi_u . H_p1 - chi_p1 . H_u
#
# H_p1=q; chi_p1=N.
# Differentiate along these directions exactly.

print("Calculating u,p1 contribution...", flush=True)

chi_u_shift = chi(u0 + e*q0, q0)
dchi_u_along_q = sp.factor(
    at_internal(sp.diff(chi_u_shift, e).subs(e, 0))
)

N0 = sp.Matrix([0, 1, 0, 0])
H_u_shift = H(u0 + e*N0, q0)
dH_u_along_N = sp.factor(
    at_internal(sp.diff(H_u_shift, e).subs(e, 0))
)

part_u = sp.factor(
    dchi_u_along_q - dH_u_along_N
)

print("chi_u . q =", dchi_u_along_q, flush=True)
print("N . H_u =", dH_u_along_N, flush=True)
print("(u,p1) =", part_u, flush=True)

# Canonical (q,P) bracket:
# chi_q . H_P - chi_P . H_q
#
# chi is independent of P, so the second term
# vanishes identically.
#
# At the constrained point H_P0=H_P1=0.
# Retain both transverse derivatives.

print("Calculating q,P contribution...", flush=True)

q_shift = sp.Matrix([a, v**2*k, z2, z3])
chi_shift = chi(u0, q_shift)

dchi_q2 = sp.factor(
    at_internal(
        sp.diff(chi_shift, z2).subs(
            {z2: 0, z3: 0}
        )
    )
)
dchi_q3 = sp.factor(
    at_internal(
        sp.diff(chi_shift, z3).subs(
            {z2: 0, z3: 0}
        )
    )
)

# From H_c, with P2=r, P3=s:
dH_dP2 = -v**5*k**2*r/2
dH_dP3 = -v**5*k**2*s/2

part_q = sp.factor(
    dchi_q2*dH_dP2
    + dchi_q3*dH_dP3
)

print("(q,P) =", part_q, flush=True)

# Internal canonical contribution:
# chi_I H_ps - chi_ps H_I.
# The toy sector has chi_ps=0.

print("Calculating internal contribution...", flush=True)

chi_internal = 2*I*k/v
H_internal = v*(ps**2/2 + I*k**2)

part_I = sp.factor(
    sp.diff(chi_internal, I)
    * sp.diff(H_internal, ps)
    -
    sp.diff(chi_internal, ps)
    * sp.diff(H_internal, I)
).subs({I: 1, ps: t})

print("(I,ps) =", part_I, flush=True)

direct = sp.factor(part_u + part_q + part_I)

print("DIRECT F_N =", direct, flush=True)

# Only now introduce the previously reconstructed
# normal-preservation expression for comparison.
proposed = (
    -n
    + 2*k*t
    - 6*k*a/v**2
    - v**3*k*cross/2
    - v**2*k*a*R/2
)

difference = sp.factor(direct - proposed)

print("PROPOSED F_N =", sp.factor(proposed), flush=True)
print("DIFFERENCE =", difference, flush=True)

if difference != 0:
    raise AssertionError(
        "Independent F_N calculation disagrees "
        "with the previously reconstructed expression"
    )

print("PASS: independent F_N family calculation")
print("GENERAL NORMAL PRESERVATION: NOT ESTABLISHED")
print("GENERAL DIRAC CLOSURE: NOT ESTABLISHED")
