"""CP098-T60: Off-surface reconstruction of E.

Compute E={chi_N,psi_T} by canonical directional
derivatives with unrestricted external momenta.

Internal sector:
  {A,B}=DA
  {A,A}=0
  {A,C}=0
where A,C depend only on internal coordinates.

No off-surface extension of the T41 result assumed.
"""

import sympy as sp

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
a = sp.symbols("a", real=True)
A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)
DA = sp.symbols("DA", real=True)
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
    w = qq + T*dot(T, qq)
    w2 = dot(w, w)
    N = w/sp.sqrt(w2)
    kap = sp.sqrt(w2)/speed**2
    return speed, T, N, kap, w2

def chiN(uu, qq):
    speed, T, N, kap, w2 = geom(uu, qq)
    return pair(p1, N)+2*A*kap/speed

def psiT(uu, qq):
    speed, T, N, kap, w2 = geom(uu, qq)
    return (
        pair(p0, T)
        + B
        - A*kap**2
        - speed**4*kap**2*dot(P,P)/(4*C)
    )

def D(expr):
    return sp.factor(sp.diff(expr, e).subs(e, 0))

print("=== CP098-T60 ===", flush=True)

speed, T, N, kap, w2 = geom(u, q)

assert sp.simplify(speed-v) == 0
assert sp.simplify(kap-k) == 0
assert sp.simplify(
    chiN(u,q)-(p1z[1]+2*A*k/v)
) == 0

print("PASS: adapted definitions", flush=True)

# Canonical bracket:
#
# {chiN,psiT} =
# chiN_q . psiT_P
# - chiN_p1 . psiT_u
# + {chiN,psiT}_internal
#
# chiN_u . psiT_p1 = 0
# chiN_P . psiT_q = 0
# No y-dependent terms.

print("Computing q/P contribution...", flush=True)

psi_P = -speed**4*kap**2*eta*P/(2*C)

term_qP = D(
    chiN(u, q+e*psi_P)
)

print("Computing p1/u contribution...", flush=True)

term_p1u = -D(
    psiT(u+e*N, q)
)

# {2*A*k/v, B}=2*k*DA/v.
# Other internal brackets vanish because
# A and C depend only on internal coordinates.

term_internal = 2*k*DA/v

external = sp.factor(term_qP+term_p1u)
E_off = sp.factor(external+term_internal)

print("q/P contribution =", term_qP, flush=True)
print("p1/u contribution =", term_p1u, flush=True)
print("Internal contribution =", term_internal,
      flush=True)
print("Unrestricted E =", E_off, flush=True)

# Independently verify the directional result
# using coordinate-gradient contraction.

q_grad = sp.Matrix([
    D(chiN(u,q+e*sp.eye(4)[:,i]))
    for i in range(4)
])

u_grad_psi = sp.Matrix([
    D(psiT(u+e*sp.eye(4)[:,i],q))
    for i in range(4)
])

gradient_external = sp.factor(
    (q_grad.T*psi_P)[0]
    -(N.T*u_grad_psi)[0]
)

assert sp.factor(
    gradient_external-external
) == 0

print("PASS: independent gradient contraction",
      flush=True)

# Impose the complete regular constraint
# family only after obtaining E_off.

r,s,x,y,n,b,c = sp.symbols(
    "r s x y n b c", real=True
)

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
        A*k**2
        +v**4*k**2*(r*r+s*s)/(4*C)
        -B,
    p0z[1]:n,
    p0z[2]:b,
    p0z[3]:c
}

E_surface = sp.factor(E_off.subs(surface))

E_T41 = -(
    4*A*C*a*k
    -4*C*DA*k*v**2
    +2*C*n*v**2
    +a*k*(r*r+s*s)*v**4
    +k*v**5*(r*x+s*y)
)/(2*C*v**3)

residual = sp.factor(E_surface-E_T41)

print("Surface E =", E_surface, flush=True)
print("T41 E =", sp.factor(E_T41), flush=True)
print("T41 residual =", residual, flush=True)

assert residual == 0

print("PASS: independent T41 cross-check",
      flush=True)

# Confirm that the unrestricted result is
# genuinely different from its surface value.

assert any(
    z in E_off.free_symbols
    for z in (Pz[0],Pz[1],p1z[0],p1z[1])
)

print("PASS: off-surface dependence retained",
      flush=True)

# Compare the off-surface E/d with its
# surface restriction. Do not assert equality.

d = 2*A/v**3
ratio_off = sp.factor(E_off/d)
ratio_surface = sp.factor(E_T41/d)

print("Off-surface E/d =", ratio_off, flush=True)
print("Surface E/d =", ratio_surface, flush=True)

assert sp.factor(
    ratio_off.subs(surface)-ratio_surface
) == 0

print("PASS: E/d surface restriction",
      flush=True)

print("OFF-SURFACE E: RECONSTRUCTED")
print("REGULAR-SURFACE T41 AGREEMENT: VERIFIED")
print("OFF-SURFACE FIRST-CLASS EXTENSION: OPEN")
print("COEFFICIENT EVOLUTION: OPEN")
print("FULL GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
