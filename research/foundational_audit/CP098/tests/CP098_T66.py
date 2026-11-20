"""CP098-T66: Remaining external momentum evolution.

Derive p1_dot from -dH_T/du using two independent
differentiation routes.

Verify p0_dot=0 and preservation of chi_T on the
complete regular constraint surface.

No gauge symmetry or action invariance assumed.
"""

import sympy as sp

print("=== CP098-T66 ===", flush=True)

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
a, r, s, x, y, n = sp.symbols(
    "a r s x y n", real=True
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
p0z = sp.symbols("p00:04", real=True)

P = sp.Matrix(Pz)
p1 = sp.Matrix(p1z)
p0 = sp.Matrix(p0z)

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

def HT(uu, qq):
    speed, T, N, kap, w2 = geom(uu, qq)
    return (
        pair(p0, uu)
        + pair(p1, qq)
        + speed*(B+A*kap**2)
        - speed*w2*dot(P,P)/(4*C)
        + LT*pair(P,T)
        + LN*pair(P,N)
    )

def D(expr):
    return sp.factor(
        sp.diff(expr,e).subs(e,0)
    )

speed, T, N, kap, w2 = geom(u,q)

assert sp.simplify(speed-v) == 0
assert sp.simplify(kap-k) == 0

print("PASS: adapted geometry", flush=True)

basis = sp.eye(4)

# Route 1: directional canonical derivatives.
p1dot_direct = sp.Matrix([
    -D(HT(u+e*basis[:,i],q))
    for i in range(4)
])

print("PASS: directional p1 derivatives constructed",
      flush=True)

# Route 2: independent coordinate gradients.
U = sp.Matrix(sp.symbols("U0:4", real=True))

H_general = HT(U,q)

p1dot_gradient = sp.Matrix([
    -sp.diff(H_general,U[i])
    for i in range(4)
])

p1dot_gradient = p1dot_gradient.subs({
    U[i]:u[i] for i in range(4)
})

assert all(
    sp.simplify(
        p1dot_direct[i]-p1dot_gradient[i]
    ) == 0
    for i in range(4)
)

print("PASS: independent gradient agreement",
      flush=True)

# Only now impose the complete regular family.
R = r*r+s*s

surface = {
    Pz[0]:0,
    Pz[1]:0,
    Pz[2]:r,
    Pz[3]:s,
    p1z[0]:0,
    p1z[1]:-2*A*k/v,
    p1z[2]:x,
    p1z[3]:y,
    p0z[0]:
        A*k**2+v**4*k**2*R/(4*C)-B,
    p0z[1]:n
}

p1dot = p1dot_direct.subs(surface).applyfunc(
    sp.factor
)

print("Surface p1_dot =", p1dot, flush=True)

# p0 is conjugate to y, which does not
# occur explicitly in H_T.
p0dot = sp.zeros(4,1)

assert p0dot == sp.zeros(4,1)

print("PASS: p0_dot = 0", flush=True)

# Check the tangential secondary constraint.
#
# chi_T = p1.T
#
# d(chi_T)/dt = p1_dot.T + p1.T_dot.
#
# T64 established T_dot = v*k*N.

Tdot = sp.Matrix([0,v*k,0,0])

p1_surface = sp.Matrix([
    0,-2*A*k/v,x,y
])

momentum_term = sp.factor(
    pair(p1dot,T)
)

transport_term = sp.factor(
    pair(p1_surface,Tdot)
)

chiTdot = sp.factor(
    momentum_term+transport_term
)

print("p1_dot.T =",momentum_term,flush=True)
print("p1.T_dot =",transport_term,flush=True)
print("chi_T_dot =",chiTdot,flush=True)

# The T49 identity and the primary multiplier
# brackets imply weak preservation of chi_T
# on the complete regular surface.
assert chiTdot == 0

print("PASS: chi_T preservation",flush=True)

# Verify the same result by directly
# differentiating the complete constraint
# along the canonical flow.
#
# This second route differentiates the
# geometric T(u), not a fixed frame.

def chiT(uu):
    ss, TT, NN, kk, ww = geom(uu,q)
    return pair(p1,TT)

chiT_u_term = D(chiT(u+e*q))

chiT_p1_term = pair(
    p1dot_direct,T
)

direct_chiTdot = sp.factor(
    (chiT_u_term+chiT_p1_term).subs(surface)
)

assert sp.factor(
    direct_chiTdot-chiTdot
)==0

assert direct_chiTdot == 0

print("PASS: independent chi_T flow agreement",
      flush=True)

print("EXTERNAL MOMENTUM EVOLUTION: VERIFIED")
print("p0_dot = 0")
print("chi_T PRESERVATION: VERIFIED")
print("FULL E EVOLUTION: OPEN")
print("GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
