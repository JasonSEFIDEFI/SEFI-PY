"""CP098-T23: Nine-parameter adapted-family bracket audit.

v,k>0; a,t,n,r,s,x,y real.
Toy internal sector A=I, C=1, B=ps**2/2.

All five constraints are checked. Canonical bracket
contributions are calculated by differentiation.
The T16 requirement uses the previously reconstructed
F_N expression, not an independent full Dirac proof.
"""

import sympy as sp

v, k = sp.symbols("v k", positive=True)
a, t, n, r, s, x, y = sp.symbols(
    "a t n r s x y", real=True
)
e, z2, z3 = sp.symbols("e z2 z3", real=True)

eta = sp.diag(-1, 1, 1, 1)

def mdot(X, Y):
    return (X.T * eta * Y)[0]

def contract(p, X):
    return (p.T * X)[0]

R = r**2 + s**2
cross = r*x + s*y

# Exact constrained point.
u0 = sp.Matrix([v, 0, 0, 0])
q0 = sp.Matrix([a, v**2*k, 0, 0])
P0 = sp.Matrix([0, 0, r, s])
p10 = sp.Matrix([0, -2*k/v, x, y])

p00 = k**2 + v**4*k**2*R/4 - t**2/2
p0 = sp.Matrix([p00, n, 0, 0])

def geometry(u, q):
    V2 = -mdot(u, u)
    T = u/sp.sqrt(V2)
    w = q + T*mdot(T, q)
    w2 = mdot(w, w)
    kap = sp.sqrt(w2)/V2
    N = w/sp.sqrt(w2)
    return V2, T, N, kap, w2

V2, T, N, kap, w2 = geometry(u0, q0)

phi_T = contract(P0, T)
phi_N = contract(P0, N)
chi_T = contract(p10, T)
chi_N = contract(p10, N) + 2*kap/v
psi_T = (
    contract(p0, T)
    + t**2/2
    - kap**2
    - w2*mdot(P0, P0)/4
)

print("=== CP098-T23 ===", flush=True)

for name, expr in [
    ("phi_T", phi_T),
    ("phi_N", phi_N),
    ("chi_T", chi_T),
    ("chi_N", chi_N),
    ("psi_T", psi_T),
]:
    value = sp.factor(expr)
    print(name, "=", value, flush=True)
    assert value == 0, name

print("PASS: all five constraints", flush=True)

# (u,p1) contribution:
# d(psi)/du along N = d(psi(u+e*N))/de at e=0.
u_e = sp.Matrix([v, e, 0, 0])
V2e, Te, Ne, ke, w2e = geometry(u_e, q0)

psi_e = (
    contract(p0, Te)
    + t**2/2
    - ke**2
    - w2e*mdot(P0, P0)/4
)

part_u = sp.factor(
    sp.diff(psi_e, e).subs(e, 0)
)
print("(u,p1) =", part_u, flush=True)

# (q,P) contribution:
# Differentiate the two transverse directions
# independently, without assuming they vanish.
q_z = sp.Matrix([a, v**2*k, z2, z3])
_, _, Nz, kz, w2z = geometry(u0, q_z)

chi_z = contract(p10, Nz) + 2*kz/v

dchi_z2 = sp.diff(chi_z, z2).subs(
    {z2: 0, z3: 0}
)
dchi_z3 = sp.diff(chi_z, z3).subs(
    {z2: 0, z3: 0}
)

# The P0,P1 derivatives vanish at the constrained
# point; chi_N is independent of P.
dpsi_dr = -w2*r/2
dpsi_ds = -w2*s/2

part_q = sp.factor(
    -dpsi_dr*dchi_z2
    -dpsi_ds*dchi_z3
)
print("(q,P) =", part_q, flush=True)

# Internal canonical bracket (I,ps).
I, ps = sp.symbols("I ps", real=True)
psi_int = ps**2/2 - I*k**2
chi_int = 2*I*k/v

part_I = sp.factor(
    sp.diff(psi_int, I)*sp.diff(chi_int, ps)
    - sp.diff(psi_int, ps)*sp.diff(chi_int, I)
).subs({I: 1, ps: t})

print("(I,ps) =", part_I, flush=True)

direct = sp.factor(part_u + part_q + part_I)

# Previously reconstructed F_N, specialized
# to this exact constrained family.
F_N = (
    -n
    + 2*k*t
    - 6*k*a/v**2
    - v**3*k*cross/2
    - v**2*k*a*R/2
)

# T16 necessary condition.
required = sp.factor(
    -4*k*a/v**3 - F_N/v
)

difference = sp.factor(direct - required)

print("Direct =", direct, flush=True)
print("T16 required =", required, flush=True)
print("Difference =", difference, flush=True)

assert sp.factor(
    part_q - v**2*k*cross/2
) == 0, "Transverse coefficient mismatch"

assert difference == 0, "T16 identity mismatch"

print("PASS: nine-parameter adapted-family identity")
print("GENERAL DIRAC CLOSURE: NOT ESTABLISHED")
