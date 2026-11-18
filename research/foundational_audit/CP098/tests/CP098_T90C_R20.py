import sympy as s

print("=== CP098-T90C-R20 ===", flush=True)

v, k = s.symbols("v k", positive=True)
A, C = s.symbols("A C", real=True)
h, j0, j1, j2, j3 = s.symbols(
    "h j0 j1 j2 j3", real=True
)
t, e = s.symbols("t e", real=True)

eta = s.diag(-1, 1, 1, 1)

def dot(x, y):
    return (x.T * eta * y)[0]

def clean(x):
    return s.simplify(x)

def geometry(u, q):
    speed = s.sqrt(-dot(u, u))
    T = u/speed
    w = q + dot(T, q)*T
    curvature = s.sqrt(dot(w, w))/speed**2
    N = w/(speed**2*curvature)
    projector = (
        s.eye(4)
        + T*(T.T*eta)
        - N*(N.T*eta)
    )
    return speed, T, curvature, N, projector

u = s.Matrix([v, 0, 0, 0])
q = s.Matrix([h, v**2*k, 0, 0])
j = s.Matrix([j0, j1, j2, j3])

speed, T, curvature, N, proj = geometry(u, q)

assert clean(speed-v) == 0
assert clean(curvature-k) == 0

R = proj*j
R2 = clean(dot(R, R))

assert clean(R2-j2**2-j3**2) == 0

# Normal variation: q -> q + e*N.
qe = q + e*N
ve, Te, ke, Ne, proje = geometry(u, qe)
Re = proje*j

Lcurv = -ve*A*ke**2
Ltors = -C*dot(Re, Re)/(ve**5*ke**2)

curv_normal = clean(
    s.diff(Lcurv, e).subs(e, 0)
)
tors_normal = clean(
    s.diff(Ltors, e).subs(e, 0)
)

# Evolve the moving frame using u'=q, q'=j.
ut = u + t*q
qt = q + t*j

vt, Tt, kt, Nt, projt = geometry(ut, qt)

Nprime = s.Matrix([
    clean(s.diff(Nt[i], t).subs(t, 0))
    for i in range(4)
])

Tprime = s.Matrix([
    clean(s.diff(Tt[i], t).subs(t, 0))
    for i in range(4)
])

assert (Tprime-v*k*N).applyfunc(clean) == s.zeros(4, 1)

normal_transverse = (proj*Nprime).applyfunc(clean)
expected_transverse = (R/(v**2*k)).applyfunc(clean)

assert (
    normal_transverse-expected_transverse
).applyfunc(clean) == s.zeros(4, 1)

P = -2*C*eta*R/(v**5*k**2)

momentum_term = clean((P.T*Nprime)[0])

expected_curvature = -2*A*k/v
expected_torsion = 2*C*R2/(v**7*k**3)

assert clean(curv_normal-expected_curvature) == 0
assert clean(tors_normal-expected_torsion) == 0
assert clean(momentum_term+expected_torsion) == 0

normal_residual = clean(
    tors_normal+momentum_term
)

assert normal_residual == 0

# Tangent variation leaves transverse acceleration
# and curvature unchanged.
q_tangent = q + e*T
vT, TT, kT, NT, projT = geometry(u, q_tangent)
RT = projT*j

Ltangent = (
    -vT*A*kT**2
    -C*dot(RT, RT)/(vT**5*kT**2)
)

tangent_gradient = clean(
    s.diff(Ltangent, e).subs(e, 0)
)

assert tangent_gradient == 0
assert clean((P.T*N)[0]) == 0

tangent_momentum_term = clean(
    (P.T*Tprime)[0]
)

assert tangent_momentum_term == 0

print("ADAPTED FRAME: arbitrary v, k, h, j0..j3")
print("R SQUARED =", R2)
print("T PRIME RESIDUAL = 0")
print("TRANSVERSE N PRIME RESIDUAL = 0")
print("TANGENT GRADIENT =", tangent_gradient)
print("TANGENT MOMENTUM TERM =", tangent_momentum_term)
print("CURVATURE NORMAL =", curv_normal)
print("TORSION NORMAL =", tors_normal)
print("MOMENTUM NORMAL TERM =", momentum_term)
print("NORMAL CANCELLATION =", normal_residual)

print()
print("R20 COMPLETE")
print("SYMBOLIC ADAPTED-FRAME IDENTITIES: VERIFIED")
print("FULL OSTROGRADSKY SYSTEM: OPEN")
print("DEFI TORSION IDENTIFICATION: UNRESOLVED")
