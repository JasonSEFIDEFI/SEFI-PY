"""CP098-T71: Off-surface p0 contribution to E.

Verify E_p0 = -(p0.N)/v from the canonical
definitions, then independently derive its
Hamiltonian evolution on the regular surface.

No full E evolution or gauge symmetry assumed.
"""

import sympy as sp

print("=== CP098-T71 ===", flush=True)

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
a, n, b, c, r, s = sp.symbols(
    "a n b c r s", real=True
)
A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)

eps = sp.symbols("eps", real=True)

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])

p0 = sp.Matrix([
    A*k**2+v**4*k**2*(r*r+s*s)/(4*C)-B,
    n, b, c
])

def dot(X, Y):
    return (X.T*eta*Y)[0]

def pair(p, X):
    return (p.T*X)[0]

def geom(uu, qq):
    speed = sp.sqrt(-dot(uu,uu))
    T = uu/speed
    w = qq+T*dot(T,qq)
    N = w/sp.sqrt(dot(w,w))
    return speed,T,N

def eq(lhs,rhs):
    return sp.cancel(sp.together(lhs-rhs)) == 0

speed,T,N = geom(u,q)

assert eq(speed,v)
assert N == sp.Matrix([0,1,0,0])

print("PASS: adapted geometry", flush=True)

# In {chi_N,psi_T}, the p0-dependent
# contribution is
#
# -sum_i N^i * d(p0.T)/du^i.
#
# Differentiate the original T(u), not a
# fixed adapted-frame expression.

basis = sp.eye(4)

canonical_term = 0

for i in range(4):
    uu = u+eps*basis[:,i]
    ss,TT,NN = geom(uu,q)

    derivative = sp.diff(
        pair(p0,TT),eps
    ).subs(eps,0)

    canonical_term -= N[i]*derivative

canonical_term = sp.factor(canonical_term)

expected = -pair(p0,N)/v

assert eq(canonical_term,expected)

print("PASS: canonical p0 contribution",
      flush=True)

print("E_p0 =",canonical_term,flush=True)

# Independently reconstruct the Hamiltonian
# evolution of E_p0=-p0.N/v.
#
# p0_dot=0 from T66.
# N_dot and v_dot come from T64.

vdot = a

Ndot = sp.Matrix([
    v*k,
    0,
    -v**3*k*r/(2*C),
    -v**3*k*s/(2*C)
])

pdot = sp.zeros(4,1)

transport_route = sp.factor(
    -(pair(pdot,N)+pair(p0,Ndot))/v
    +pair(p0,N)*vdot/v**2
)

# Independent route: perturb u and q by
# their canonical velocities and differentiate
# the full geometric expression -p0.N/v.

LT,LN = sp.symbols(
    "lambda_T lambda_N",real=True
)

udot = q

P = sp.Matrix([0,0,r,s])

qdot = sp.Matrix([
    LT,
    LN,
    -v**5*k**2*r/(2*C),
    -v**5*k**2*s/(2*C)
])

ue = u+eps*udot
qe = q+eps*qdot

ve,Te,Ne = geom(ue,qe)

full_expression = -pair(p0,Ne)/ve

directional_route = sp.factor(
    sp.diff(full_expression,eps).subs(eps,0)
)

assert eq(
    transport_route,directional_route
)

print("PASS: independent geometric evolution",
      flush=True)

R = r*r+s*s

expected_dot = (
    -k*p0[0]
    +v**2*k*(b*r+c*s)/(2*C)
    +a*n/v**2
)

assert eq(
    transport_route,expected_dot
)

print("PASS: transverse p0 terms retained",
      flush=True)

print(
    "E_p0_dot =",
    sp.factor(transport_route),
    flush=True
)

# The result must be independent of the
# unrestricted tangential multiplier and
# of the normal multiplier.

assert LT not in transport_route.free_symbols
assert LN not in transport_route.free_symbols

print("PASS: multiplier independence",
      flush=True)

# Verify the explicit transverse coefficients.
assert eq(
    sp.diff(transport_route,b),
    v**2*k*r/(2*C)
)

assert eq(
    sp.diff(transport_route,c),
    v**2*k*s/(2*C)
)

print("PASS: b,c derivative coefficients",
      flush=True)

print("P0 CONTRIBUTION TO E: VERIFIED")
print("P0 CONTRIBUTION TO E_DOT: VERIFIED")
print("FULL E EVOLUTION: OPEN")
print("JACOBI CROSS-CHECK: OPEN")
print("GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
