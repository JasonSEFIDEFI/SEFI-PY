"""CP098-T85A: General geometric reconstruction of E.

Construct {chi_N,psi_T} without imposing an
adapted frame during its geometric derivation.

Then compare with the independently established
T70/T60 adapted expression.

Does NOT calculate {E,H_T}.
"""

import sympy as sp

print("=== CP098-T85A ===", flush=True)

v, k = sp.symbols("v k", positive=True)
a = sp.symbols("a", real=True)
A, C = sp.symbols("A C", real=True, nonzero=True)
B, DA = sp.symbols("B DA", real=True)

eta = sp.diag(-1, 1, 1, 1)
eps = sp.symbols("eps", real=True)

u = sp.Matrix(sp.symbols("u0:4", real=True))
q = sp.Matrix(sp.symbols("q0:4", real=True))
P = sp.Matrix(sp.symbols("P0:4", real=True))
p1 = sp.Matrix(sp.symbols("j0:4", real=True))
p0 = sp.Matrix(sp.symbols("z0:4", real=True))

def dot(X, Y):
    return (X.T*eta*Y)[0]

def pair(p, X):
    return (p.T*X)[0]

def eq(f, g):
    return sp.cancel(sp.together(f-g)) == 0

def geom(uu, qq):
    speed = sp.sqrt(-dot(uu, uu))
    T = uu/speed
    w = qq+T*dot(T, qq)
    length = sp.sqrt(dot(w,w))
    N = w/length
    kap = length/speed**2
    return speed, T, N, kap, length

speed, T, N, kap, length = geom(u,q)

# The derivative of N with respect to q
# follows from differentiating its normalized
# projected-vector definition.
#
# dN^mu/dq^i =
# (delta^mu_i + T^mu T_i - N^mu N_i)/length.
#
# All canonical momenta are covectors.

chi_q = sp.Matrix([
    (
        p1[i]
        +pair(p1,T)*(eta*T)[i]
        -pair(p1,N)*(eta*N)[i]
    )/length
    +2*A*(eta*N)[i]/speed**3
    for i in range(4)
])

psi_P = (
    -speed**4*kap**2*(eta*P)/(2*C)
)

# Calculate N^i d(psi_T)/du^i using a
# directional derivative with N held fixed
# as the direction at the base point.

shifted_u = u+eps*N
sv, Tv, Nv, kv, lv = geom(shifted_u,q)

psi_shifted = (
    pair(p0,Tv)
    +B
    -A*kv**2
    -sv**4*kv**2*dot(P,P)/(4*C)
)

normal_psi_u = sp.diff(
    psi_shifted,eps
).subs(eps,0)

E_general = (
    (chi_q.T*psi_P)[0]
    -normal_psi_u
    +2*kap*DA/speed
)

print("PASS: unrestricted geometric bracket assembled",
      flush=True)

# Only now evaluate at the adapted point.

adapted = {
    u[0]:v, u[1]:0, u[2]:0, u[3]:0,
    q[0]:a, q[1]:v**2*k, q[2]:0, q[3]:0
}

E_adapted = sp.factor(
    E_general.subs(adapted)
)

print("PASS: adapted evaluation completed",
      flush=True)

P0,P1,P2,P3 = P
j0,j1,j2,j3 = p1
z0,z1,z2,z3 = p0

P_squared = (
    -P0**2+P1**2+P2**2+P3**2
)

E_T70 = -(
    4*A*C*a*k
    +2*A*P1*k**2*v**4
    -4*C*DA*k*v**2
    +2*C*z1*v**2
    +a*k*v**4*P_squared
    +k*v**5*(P2*j2+P3*j3)
)/(2*C*v**3)

residual = sp.factor(E_adapted-E_T70)

print("T70_RESIDUAL =",residual,flush=True)

assert residual == 0, (
    "STOP: general bracket disagrees with T70"
)

print("PASS: unrestricted bracket reproduces T70",
      flush=True)

# Verify that the regular surface has not
# been imposed during the bracket calculation.

assert P0 in E_adapted.free_symbols
assert P1 in E_adapted.free_symbols

print("PASS: off-surface momenta retained",
      flush=True)

r,s,x,y,n = sp.symbols(
    "r s x y n", real=True
)

surface = {
    P0:0, P1:0, P2:r, P3:s,
    j0:0, j1:-2*A*k/v,
    j2:x, j3:y,
    z1:n
}

E_surface = sp.factor(
    E_adapted.subs(surface)
)

expected = (
    -2*A*a*k/v**3
    +2*DA*k/v
    -n/v
    -a*k*v*(r**2+s**2)/(2*C)
    -k*v**2*(r*x+s*y)/(2*C)
)

assert eq(E_surface,expected), (
    "STOP: regular-surface E mismatch"
)

print("PASS: T74 surface E recovered",
      flush=True)

print("GENERAL GEOMETRIC E: RECONSTRUCTED")
print("ADAPTED T70 AGREEMENT: EXACT")
print("REGULAR T74 E AGREEMENT: EXACT")
print("FULL CANONICAL E_DOT: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T85A COMPLETE")
