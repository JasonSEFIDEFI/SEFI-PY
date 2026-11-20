"""CP098-T72: p1-dependent contribution to E.

1. Derive E_p1 from canonical directional derivatives.
2. Verify its covariant geometric expression.
3. Compute its regular-surface evolution by two routes.
4. Retain moving-frame and internal C evolution.

This is not yet the full E evolution.
"""

import sympy as sp

print("=== CP098-T72 ===", flush=True)

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
a, r, s, x, y, n, b, c = sp.symbols(
    "a r s x y n b c", real=True
)
A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)
DC = sp.symbols("DC", real=True)
LT, LN = sp.symbols("LT LN", real=True)
e = sp.symbols("e", real=True)

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])
P = sp.Matrix([0, 0, r, s])
p1 = sp.Matrix([0, -2*A*k/v, x, y])

R = r*r+s*s
S = r*x+s*y

def dot(X, Y):
    return (X.T*eta*Y)[0]

def pair(p, X):
    return (p.T*X)[0]

def eq(X, Y):
    return sp.cancel(sp.together(X-Y)) == 0

def geom(uu, qq):
    speed = sp.sqrt(-dot(uu,uu))
    T = uu/speed
    w = qq+T*dot(T,qq)
    wn = sp.sqrt(dot(w,w))
    N = w/wn
    kap = wn/speed**2
    return speed,T,N,kap

def Wfun(pp, PP, TT, NN):
    return (
        pair(pp,eta*PP)
        +pair(pp,TT)*pair(PP,TT)
        -pair(pp,NN)*pair(PP,NN)
    )

speed,T,N,kap = geom(u,q)

assert eq(speed,v)
assert eq(kap,k)
assert eq(Wfun(p1,P,T,N),S)

print("PASS: regular adapted geometry",flush=True)

# First verify the canonical origin of E_p1
# without imposing momentum constraints.

Pz = sp.Matrix(sp.symbols("P0:4",real=True))
jz = sp.Matrix(sp.symbols("j0:4",real=True))

basis = sp.eye(4)

canonical = 0

for i in range(4):
    qe = q+e*basis[:,i]
    se,Te,Ne,ke = geom(u,qe)

    derivative = sp.diff(
        pair(jz,Ne),e
    ).subs(e,0)

    dpsi_dP = (
        -v**4*k**2*(eta*Pz)[i]/(2*C)
    )

    canonical += derivative*dpsi_dP

geometric = (
    -v**2*k*Wfun(jz,Pz,T,N)/(2*C)
)

assert eq(canonical,geometric)

print("PASS: off-surface canonical identity",
      flush=True)

# Now use the complete regular constraint
# family and its verified canonical flow.

udot = q

qdot = sp.Matrix([
    LT,
    LN,
    -v**5*k**2*r/(2*C),
    -v**5*k**2*s/(2*C)
])

Pdot = sp.Matrix([
    0,
    v**3*k*R/(2*C),
    -x-LN*r/(v**2*k),
    -y-LN*s/(v**2*k)
])

p1dot = sp.Matrix([
    2*A*k**2,
    2*A*a*k/v**2-n-a*k*v**2*R/(2*C),
    a*LN*r/(k*v**3)-LT*r/v-b,
    a*LN*s/(k*v**3)-LT*s/v-c
])

Tdot = sp.Matrix([0,v*k,0,0])

Ndot = sp.Matrix([
    v*k,
    0,
    -v**3*k*r/(2*C),
    -v**3*k*s/(2*C)
])

vdot = a
kdot = LN/v**2-3*a*k/v
Cdot = v*DC

# Route 1: complete product-rule transport.

def product_dot(p,pdot,X,Xdot):
    return pair(pdot,X)+pair(p,Xdot)

Wdot = (
    product_dot(p1,p1dot,eta*P,eta*Pdot)
    +product_dot(p1,p1dot,T,Tdot)*pair(P,T)
    +pair(p1,T)*product_dot(P,Pdot,T,Tdot)
    -product_dot(p1,p1dot,N,Ndot)*pair(P,N)
    -pair(p1,N)*product_dot(P,Pdot,N,Ndot)
)

Wdot = sp.factor(Wdot)

prefactor = -v**2*k/(2*C)

prefactor_dot = (
    -v*a*k/C
    -v**2*kdot/(2*C)
    +v**2*k*Cdot/(2*C**2)
)

E_dot_transport = sp.factor(
    prefactor_dot*S+prefactor*Wdot
)

print("PASS: product-rule evolution constructed",
      flush=True)

# Route 2: evolve the geometric expression
# itself before evaluating the derivative.
#
# The frame is recomputed from u(e),q(e).
# All momenta and C evolve simultaneously.

ue = u+e*udot
qe = q+e*qdot
Pe = P+e*Pdot
p1e = p1+e*p1dot
Ce = C+e*Cdot

ve,Te,Ne,ke = geom(ue,qe)

We = Wfun(p1e,Pe,Te,Ne)

# Compare the independent derivative of W.
Wdot_direct = sp.diff(We,e).subs(e,0)

assert eq(Wdot_direct,Wdot)

print("PASS: independent moving-frame W evolution",
      flush=True)

# Differentiate the entire geometric expression.
E_e = -ve**2*ke*We/(2*Ce)

E_dot_direct = sp.diff(E_e,e).subs(e,0)

assert eq(E_dot_direct,E_dot_transport)

print("PASS: independent E_p1 evolution",
      flush=True)

# Ensure internal C evolution is retained.
expected_DC_coefficient = (
    v**3*k*S/(2*C**2)
)

assert eq(
    sp.diff(E_dot_transport,DC),
    expected_DC_coefficient
)

print("PASS: C-gradient contribution retained",
      flush=True)

# The full regular family includes b and c.
# They enter through p1dot and must survive.
assert eq(
    sp.diff(E_dot_transport,b),
    v**2*k*r/(2*C)
)

assert eq(
    sp.diff(E_dot_transport,c),
    v**2*k*s/(2*C)
)

print("PASS: transverse p0 contributions retained",
      flush=True)

print("W_dot =",Wdot,flush=True)
print(
    "E_p1_dot =",
    sp.factor(E_dot_transport),
    flush=True
)

print("P1 CONTRIBUTION TO E: VERIFIED")
print("P1 CONTRIBUTION TO E_DOT: VERIFIED")
print("FULL E EVOLUTION: OPEN")
print("JACOBI CROSS-CHECK: OPEN")
print("GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
