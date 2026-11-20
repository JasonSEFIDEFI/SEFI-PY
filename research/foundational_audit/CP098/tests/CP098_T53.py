"""CP098-T53: Primary first-class canonical flow.

Construct the Hamiltonian vector field of phi_T
and differentiate the constraints and Hc along it.

All external momenta remain unrestricted until
the final surface evaluation.

No full gauge generator is assumed.
"""

import sympy as sp

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
a = sp.symbols("a", real=True)
A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)
e = sp.symbols("e", real=True)

Pz = sp.symbols("P0:4", real=True)
p1z = sp.symbols("p10:14", real=True)
p0z = sp.symbols("p00:04", real=True)

P = sp.Matrix(Pz)
p1 = sp.Matrix(p1z)
p0 = sp.Matrix(p0z)

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])

def dot(X, Y):
    return (X.T*eta*Y)[0]

def pair(p, X):
    return (p.T*X)[0]

def geom(uu, qq):
    speed = sp.sqrt(-dot(uu, uu))
    T = uu/speed
    w = qq+T*dot(T, qq)
    w2 = dot(w, w)
    N = w/sp.sqrt(w2)
    kap = sp.sqrt(w2)/speed**2
    return speed, T, N, kap, w2

def constraints(uu, qq, pp1):
    speed, T, N, kap, w2 = geom(uu, qq)
    return sp.Matrix([
        pair(P, T),
        pair(P, N),
        pair(pp1, T),
        pair(pp1, N)+2*A*kap/speed,
        pair(p0, T)+B-A*kap**2
        -speed**4*kap**2*dot(P,P)/(4*C)
    ])

def H(uu, qq, pp1):
    speed, T, N, kap, w2 = geom(uu, qq)
    return (
        pair(p0, uu)
        +pair(pp1, qq)
        +speed*(B+A*kap**2)
        -speed*w2*dot(P,P)/(4*C)
    )

def D(expr):
    return sp.factor(sp.diff(expr,e).subs(e,0))

print("=== CP098-T53 ===", flush=True)

speed, T, N, kap, w2 = geom(u,q)
phiT = pair(P,T)

# For G=phi_T=P_mu T^mu:
#
# delta q^mu = {q^mu,G} = T^mu
# delta p1_mu = {p1_mu,G} = -dG/du^mu
# delta u = delta P = delta y = delta p0 = 0
#
# dG/du^mu = (P_mu+G*T_mu)/v.

delta_q = T
delta_p1 = -(P+phiT*eta*T)/speed

expected_delta_p1 = sp.Matrix([
    0, -Pz[1]/v, -Pz[2]/v, -Pz[3]/v
])

assert all(
    sp.simplify(z)==0
    for z in delta_p1-expected_delta_p1
)

print("delta q =", delta_q.T, flush=True)
print("delta p1 =", delta_p1.T, flush=True)
print("PASS: canonical vector field", flush=True)

# Independently differentiate all constraints along
# the canonical vector field.

old_constraints = constraints(u,q,p1)
shifted_constraints = constraints(
    u,
    q+e*delta_q,
    p1+e*delta_p1
)

variations = sp.Matrix([
    D(shifted_constraints[i])
    for i in range(5)
])

names = ["phi_T","phi_N","chi_T","chi_N","psi_T"]

print("Constraint variations:")
for name, val in zip(names,variations):
    print(name, ":", val, flush=True)

expected = sp.Matrix([
    0, 0, 0, -Pz[1]/v, 0
])

assert all(
    sp.factor(variations[i]-expected[i])==0
    for i in range(5)
)

print("PASS: exact adapted-frame variations",
      flush=True)

# Check the Hamiltonian variation directly,
# without using T46's bracket as an input.

delta_H = D(H(
    u,
    q+e*delta_q,
    p1+e*delta_p1
))

expected_H = p1z[0]-v*k*Pz[1]

print("delta Hc =", delta_H, flush=True)

assert sp.factor(delta_H-expected_H)==0

print("PASS: independent Hamiltonian-flow check",
      flush=True)

# Tangency test: impose only the primary
# constraints P0=P1=0.
primary_surface = {
    Pz[0]:0,
    Pz[1]:0
}

primary_variations = variations.subs(primary_surface)

assert primary_variations == sp.zeros(5,1)

print(
    "PASS: constraint variations vanish "
    "on primary surface",
    flush=True
)

# Hamiltonian variation requires the secondary
# condition chi_T=0 as well.
assert sp.factor(
    delta_H.subs(primary_surface)-p1z[0]
)==0

assert sp.factor(
    delta_H.subs({
        Pz[0]:0,
        Pz[1]:0,
        p1z[0]:0
    })
)==0

print(
    "PASS: Hc variation vanishes on "
    "primary-secondary surface",
    flush=True
)

# Confirm that the flow is not identically zero:
# transverse p1 momenta change if transverse
# P components are nonzero.
assert delta_p1[2] != 0
assert delta_p1[3] != 0

print("PASS: nontrivial canonical flow", flush=True)

print("PRIMARY FIRST-CLASS FLOW: VERIFIED")
print("REGULAR CONSTRAINT-SURFACE TANGENCY: VERIFIED")
print("FULL TIME-DEPENDENT GENERATOR: OPEN")
print("CANONICAL ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
