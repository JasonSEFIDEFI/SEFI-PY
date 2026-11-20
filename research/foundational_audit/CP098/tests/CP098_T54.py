"""CP098-T54: Exact finite phi_T Hamiltonian flow.

Verify the finite canonical transformation generated
by phi_T, its composition law, and its action on
the five regular constraints and Hc.

No time-dependent gauge symmetry is assumed.
"""

import sympy as sp

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
a, eps, eps2 = sp.symbols(
    "a eps eps2", real=True
)
A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)

P0, P1, P2, P3 = sp.symbols(
    "P0 P1 P2 P3", real=True
)
p10, p11, p12, p13 = sp.symbols(
    "p10 p11 p12 p13", real=True
)
p00, p01, p02, p03 = sp.symbols(
    "p00 p01 p02 p03", real=True
)

P = sp.Matrix([P0, P1, P2, P3])
p1 = sp.Matrix([p10, p11, p12, p13])
p0 = sp.Matrix([p00, p01, p02, p03])

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])

def dot(X, Y):
    return (X.T*eta*Y)[0]

def pair(p, X):
    return (p.T*X)[0]

def geometry(uu, qq):
    speed = sp.sqrt(-dot(uu, uu))
    T = uu/speed
    w = qq+T*dot(T, qq)
    w2 = dot(w, w)
    N = w/sp.sqrt(w2)
    kap = sp.sqrt(w2)/speed**2
    return speed, T, N, kap, w2

def constraints(uu, qq, pp1):
    speed, T, N, kap, w2 = geometry(uu, qq)
    return sp.Matrix([
        pair(P, T),
        pair(P, N),
        pair(pp1, T),
        pair(pp1, N)+2*A*kap/speed,
        pair(p0, T)+B-A*kap**2
        -speed**4*kap**2*dot(P,P)/(4*C)
    ])

def H(uu, qq, pp1):
    speed, T, N, kap, w2 = geometry(uu, qq)
    return (
        pair(p0, uu)
        + pair(pp1, qq)
        + speed*(B+A*kap**2)
        - speed*w2*dot(P,P)/(4*C)
    )

def zero(expr):
    return sp.factor(expr) == 0

print("=== CP098-T54 ===", flush=True)

speed, T, N, kap, w2 = geometry(u,q)
phi_T = pair(P,T)

# The Hamiltonian vector field from T53:
#
# dq/ds = T
# dp1/ds = -(P+phi_T eta T)/v
#
# u, P, p0 and the internal variables are constant.
# Consequently T and the vector field coefficients
# remain constant along this flow.

K = (P+phi_T*eta*T)/v

def flow(qq, pp1, parameter):
    return (
        qq+parameter*T,
        pp1-parameter*K
    )

q_new, p1_new = flow(q,p1,eps)

print("Finite q =", q_new.T, flush=True)
print("Finite p1 =", p1_new.T, flush=True)

# Group composition:
q_twice, p1_twice = flow(q_new,p1_new,eps2)
q_sum, p1_sum = flow(q,p1,eps+eps2)

assert all(
    zero(z) for z in q_twice-q_sum
)
assert all(
    zero(z) for z in p1_twice-p1_sum
)

print("PASS: exact additive composition law",
      flush=True)

# Geometry should remain unchanged except for
# the tangential component of q.
new_speed, new_T, new_N, new_k, new_w2 = (
    geometry(u,q_new)
)

assert zero(new_speed-speed)
assert all(zero(z) for z in new_T-T)
assert all(zero(z) for z in new_N-N)
assert zero(new_k-k)
assert zero(new_w2-w2)

print("PASS: regular geometry preserved",
      flush=True)

old_constraints = constraints(u,q,p1)
new_constraints = constraints(u,q_new,p1_new)

change = (new_constraints-old_constraints).applyfunc(
    sp.factor
)

print("Finite constraint changes:")
for label, value in zip(
    ["phi_T","phi_N","chi_T","chi_N","psi_T"],
    change
):
    print(label, ":", value, flush=True)

expected = sp.Matrix([
    0,
    0,
    0,
    -eps*P1/v,
    0
])

assert all(
    zero(change[i]-expected[i])
    for i in range(5)
)

print("PASS: exact finite constraint transformation",
      flush=True)

# Hamiltonian transformation.
H_change = sp.factor(
    H(u,q_new,p1_new)-H(u,q,p1)
)

expected_H = eps*(p10-v*k*P1)

print("Finite Hamiltonian change =", H_change,
      flush=True)

assert zero(H_change-expected_H)

print("PASS: exact finite Hamiltonian identity",
      flush=True)

# On the full constraint surface:
# P0=P1=p10=0, p11=-2*A*k/v,
# and psi_T fixes p00.
#
# Only P1 and p10 are needed to make the
# finite changes vanish.
surface = {
    P0:0,
    P1:0,
    p10:0
}

assert all(
    zero(z.subs(surface)) for z in change
)

assert zero(H_change.subs(surface))

print("PASS: complete constraint surface invariant",
      flush=True)
print("PASS: Hc invariant on constraint surface",
      flush=True)

# The flow is nontrivial even on the surface
# when transverse P is nonzero.
assert sp.factor(
    (p1_new[2]-p1[2]).subs(surface)
    +eps*P2/v
)==0

print("PASS: nontrivial finite flow", flush=True)

print("FINITE PRIMARY HAMILTONIAN FLOW: VERIFIED")
print("TIME-DEPENDENT GAUGE GENERATOR: OPEN")
print("CANONICAL ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
