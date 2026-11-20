"""CP098-T70: Independent reconstruction of E.

Compute {chi_N,psi_T} at an adapted geometric
point using directional canonical derivatives.

All external momenta remain unrestricted.
Internal bracket is retained as DA={A,B}.

This verifies E itself, not its time derivative.
"""

import sympy as sp

print("=== CP098-T70 ===", flush=True)

v, k = sp.symbols("v k", positive=True)
a = sp.symbols("a", real=True)
A, C = sp.symbols("A C", real=True, nonzero=True)
B, DA = sp.symbols("B DA", real=True)

eta = sp.diag(-1, 1, 1, 1)
eps = sp.symbols("eps", real=True)

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])

P = sp.Matrix(sp.symbols("P0:4", real=True))
p1 = sp.Matrix(sp.symbols("j0:4", real=True))
p0 = sp.Matrix(sp.symbols("z0:4", real=True))

basis = sp.eye(4)

def dot(X, Y):
    return (X.T*eta*Y)[0]

def pair(p, X):
    return (p.T*X)[0]

def geom(uu, qq):
    speed = sp.sqrt(-dot(uu, uu))
    T = uu/speed
    w = qq+T*dot(T, qq)
    kap = sp.sqrt(dot(w,w))/speed**2
    N = w/sp.sqrt(dot(w,w))
    return speed, T, N, kap

def chiN(uu, qq):
    speed, T, N, kap = geom(uu, qq)
    return pair(p1,N)+2*A*kap/speed

def psiT(uu, qq):
    speed, T, N, kap = geom(uu, qq)
    return (
        pair(p0,T)
        +B
        -A*kap**2
        -speed**4*kap**2*dot(P,P)/(4*C)
    )

def directional(func, which, i):
    if which == "u":
        expression = func(u+eps*basis[:,i],q)
    else:
        expression = func(u,q+eps*basis[:,i])

    return sp.diff(expression,eps).subs(eps,0)

speed, T, N, kap = geom(u,q)

assert sp.simplify(speed-v)==0
assert sp.simplify(kap-k)==0

print("PASS: adapted evaluation geometry",
      flush=True)

# Canonical pairs:
# (u,p1), (q,P), (y,p0), (I,pI).
#
# chiN is independent of P and pI.
# psiT is independent of p1.
#
# Therefore:
#
# E_external =
# sum_i [
#   (d chiN/dq_i)(d psiT/dP_i)
#   -(d chiN/dp1_i)(d psiT/du_i)
# ]
#
# d chiN/dp1_i = N^i.
#
# E_internal = 2*k*DA/v.

E_external = 0

for i in range(4):
    dq_chi = directional(chiN,"q",i)

    dP_psi = (
        -v**4*k**2*(eta*P)[i]/(2*C)
    )

    du_psi = directional(psiT,"u",i)

    E_external += (
        dq_chi*dP_psi
        -N[i]*du_psi
    )

print("PASS: external directional bracket built",
      flush=True)

E_internal = 2*k*DA/v

E_direct = sp.factor(E_external+E_internal)

print("PASS: internal DA contribution included",
      flush=True)

P0,P1,P2,P3 = P
j0,j1,j2,j3 = p1
z0,z1,z2,z3 = p0

P_squared = -P0**2+P1**2+P2**2+P3**2

# Independent T60 expression.
E_T60 = -(
    4*A*C*a*k
    +2*A*P1*k**2*v**4
    -4*C*DA*k*v**2
    +2*C*z1*v**2
    +a*k*v**4*P_squared
    +k*v**5*(P2*j2+P3*j3)
)/(2*C*v**3)

assert sp.factor(E_direct-E_T60)==0

print("PASS: exact agreement with T60",
      flush=True)

# No regular momentum constraints have been
# imposed in the preceding comparison.
assert P0 in E_direct.free_symbols
assert P1 in E_direct.free_symbols

print("PASS: off-surface momenta retained",
      flush=True)

# Only now impose the regular constraint family.

r,s,x,y,n = sp.symbols(
    "r s x y n", real=True
)

surface = {
    P0:0,
    P1:0,
    P2:r,
    P3:s,
    j0:0,
    j1:-2*A*k/v,
    j2:x,
    j3:y,
    z1:n
}

E_surface = sp.factor(E_direct.subs(surface))

expected_surface = -(
    4*A*C*a*k
    -4*C*DA*k*v**2
    +2*C*n*v**2
    +a*k*v**4*(r**2+s**2)
    +k*v**5*(r*x+s*y)
)/(2*C*v**3)

assert sp.factor(
    E_surface-expected_surface
)==0

print("PASS: regular-surface E recovered",
      flush=True)

print("E_surface =",E_surface,flush=True)

print("CANONICAL E RECONSTRUCTION: VERIFIED")
print("T60 EXPRESSION: INDEPENDENTLY REPRODUCED")
print("REGULAR SURFACE: APPLIED AFTER BRACKET")
print("FULL E EVOLUTION: OPEN")
print("JACOBI CROSS-CHECK: OPEN")
print("GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
