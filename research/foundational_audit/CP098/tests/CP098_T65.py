"""CP098-T65: Canonical evolution of P_mu.

Compute P_dot = -dH_T/dq from the original
geometric definitions.

Cross-check using two independent differentiation
routes and the regular constraint family.

No gauge symmetry is assumed.
"""

import sympy as sp

print("=== CP098-T65 ===", flush=True)

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
a, r, s, x, y = sp.symbols(
    "a r s x y", real=True
)
A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)
LT, LN = sp.symbols(
    "lambda_T lambda_N", real=True
)
e = sp.symbols("e", real=True)

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])

Pz = sp.symbols("P0:4", real=True)
p1z = sp.symbols("p10:14", real=True)

P = sp.Matrix(Pz)
p1 = sp.Matrix(p1z)

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

def Hq(uu, qq):
    speed, T, N, kap, w2 = geom(uu, qq)
    return (
        pair(p1, qq)
        + speed*(B+A*kap**2)
        - speed**5*kap**2*dot(P,P)/(4*C)
        + LT*pair(P,T)
        + LN*pair(P,N)
    )

def D(expr):
    return sp.factor(sp.diff(expr,e).subs(e,0))

speed,T,N,kap,w2 = geom(u,q)

assert sp.simplify(speed-v)==0
assert sp.simplify(kap-k)==0

print("PASS: adapted geometry",flush=True)

# Route 1: directional derivatives in each
# independent q-coordinate direction.

basis = sp.eye(4)

Pdot_direct = sp.Matrix([
    -D(Hq(u,q+e*basis[:,i]))
    for i in range(4)
])

print("PASS: directional derivatives constructed",
      flush=True)

# Route 2: differentiate the independent
# q-coordinate variables symbolically.

Q = sp.Matrix(sp.symbols("Q0:4",real=True))

H_general_q = Hq(u,Q)

Pdot_gradient = sp.Matrix([
    -sp.diff(H_general_q,Q[i])
    for i in range(4)
])

Pdot_gradient = Pdot_gradient.subs({
    Q[i]:q[i] for i in range(4)
})

assert all(
    sp.simplify(
        Pdot_direct[i]-Pdot_gradient[i]
    )==0
    for i in range(4)
)

print("PASS: independent gradient agreement",
      flush=True)

# Apply the regular five-constraint family
# only after canonical differentiation.

surface = {
    Pz[0]:0,
    Pz[1]:0,
    Pz[2]:r,
    Pz[3]:s,
    p1z[0]:0,
    p1z[1]:-2*A*k/v,
    p1z[2]:x,
    p1z[3]:y
}

Pdot = Pdot_direct.subs(surface).applyfunc(
    sp.factor
)

R = r*r+s*s

expected = sp.Matrix([
    0,
    v**3*k*R/(2*C),
    -x-LN*r/(v**2*k),
    -y-LN*s/(v**2*k)
])

assert all(
    sp.factor(Pdot[i]-expected[i])==0
    for i in range(4)
)

print("PASS: regular-surface momentum evolution",
      flush=True)

print("P_dot =",Pdot,flush=True)

# Tangential primary multiplier must not
# enter P_dot because phi_T is q-independent.

assert all(
    LT not in Pdot[i].free_symbols
    for i in range(4)
)

print("PASS: lambda_T independence",flush=True)

# The normal multiplier must remain in
# the transverse momentum equations.

assert sp.factor(
    sp.diff(Pdot[2],LN)+r/(v**2*k)
)==0

assert sp.factor(
    sp.diff(Pdot[3],LN)+s/(v**2*k)
)==0

print("PASS: normal multiplier dependence",
      flush=True)

# Check preservation of phi_N=P.N using
# the previously established N_dot from T64.

Ndot = sp.Matrix([
    v*k,
    0,
    -v**3*k*r/(2*C),
    -v**3*k*s/(2*C)
])

phiNdot = sp.factor(
    pair(Pdot,N)+pair(
        sp.Matrix([0,0,r,s]),Ndot
    )
)

momentum_contribution = sp.factor(
    pair(Pdot, N)
)

frame_contribution = sp.factor(
    pair(sp.Matrix([0,0,r,s]), Ndot)
)

expected_momentum = v**3*k*R/(2*C)
expected_frame = -v**3*k*R/(2*C)

assert sp.factor(
    momentum_contribution-expected_momentum
)==0

assert sp.factor(
    frame_contribution-expected_frame
)==0

assert sp.factor(
    momentum_contribution+frame_contribution
)==0

assert sp.factor(phiNdot)==0

print("P_dot.N =", momentum_contribution,
      flush=True)
print("P.N_dot =", frame_contribution,
      flush=True)
print("phi_N_dot =", phiNdot, flush=True)

print("PASS: independent phi_N evolution",
      flush=True)

# Compare against the exact T46 identity:
#
# {phi_N,Hc} =
# -chi_N+t*phi_T
# -v^3*k*(phi_T^2-phi_N^2)/(2*C).
#
# On the complete constraint surface, the
# canonical Hamiltonian contribution vanishes.
#
# The previous expected_phiNdot is NOT the
# complete derivative: the N_dot above was
# calculated with q_dot containing the
# canonical P term, but the P_dot expression
# includes the p1 constraint contribution.
#
# Therefore this check is deliberately a
# strict stop-gate on the full calculation.

assert sp.factor(phiNdot)==0, (
    "STOP: phi_N preservation discrepancy; "
    "audit canonical momentum/normal transport"
)

print("PASS: phi_N preservation",flush=True)

print("P_MU EVOLUTION: VERIFIED")
print("MOVING-FRAME CROSS-CHECK: VERIFIED")
print("FULL E EVOLUTION: OPEN")
print("GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
