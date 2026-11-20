import sympy as s

print("=== CP098-T90C-R22 ===", flush=True)

eta = s.diag(-1, 1, 1, 1)
t, e = s.symbols("t e", real=True)

def dot(x, y):
    return (x.T * eta * y)[0]

def clean(x):
    return s.simplify(x)

def geom(u, q):
    v = s.sqrt(-dot(u, u))
    T = u/v
    w = q + dot(T, q)*T
    k = s.sqrt(dot(w, w))/v**2
    N = w/(v**2*k)
    proj = s.eye(4) + T*(T.T*eta) - N*(N.T*eta)
    return v, T, k, N, proj

def L(u, q, j, A, C):
    v, T, k, N, proj = geom(u, q)
    R = proj*j
    return -v*A*k**2 - C*dot(R, R)/(v**5*k**2)

def P(u, q, j, C):
    v, T, k, N, proj = geom(u, q)
    return -2*C*eta*proj*j/(v**5*k**2)

u = s.Matrix([5, 3, 0, 0])
q = s.Matrix([1, 2, 2, 1])
j = s.Matrix([4, 1, 3, 6])
a = s.Matrix([2, -1, 4, 3])
b = s.Matrix([1, 3, -2, 5])

A0, Adot = s.Integer(3), s.Integer(2)
C0, Cdot, Cddot = s.Integer(2), s.Integer(1), s.Integer(3)

ut = u + q*t + j*t**2/2
qt = q + j*t + a*t**2/2
jt = j + a*t + b*t**2/2
At = A0 + Adot*t
Ct = C0 + Cdot*t + Cddot*t**2/2

v, T, k, N, proj = geom(u, q)
print("v =", v, "kappa =", k, flush=True)

# First contribution: partial L / partial u at t=0.
print("STAGE 1: partial L / partial u", flush=True)
gu = []
for i in range(4):
    basis = s.eye(4)[:, i]
    expr = L(u + e*basis, q, j, A0, C0)
    gu.append(clean(s.diff(expr, e).subs(e, 0)))
gu = s.Matrix(gu)
print("STAGE 1 COMPLETE", flush=True)

# Second contribution: total derivative of partial L / partial q.
# Differentiate in t before substituting t=0.
print("STAGE 2: d/dlambda(partial L / partial q)", flush=True)
gqdot = []
for i in range(4):
    basis = s.eye(4)[:, i]
    expr = L(ut, qt + e*basis, jt, At, Ct)
    mixed = s.diff(s.diff(expr, e), t)
    gqdot.append(clean(mixed.subs({e: 0, t: 0})))
gqdot = s.Matrix(gqdot)
print("STAGE 2 COMPLETE", flush=True)

# Third contribution: second total derivative of P.
print("STAGE 3: d2P/dlambda2", flush=True)
Pt = P(ut, qt, jt, Ct)
Pddot = s.Matrix([
    clean(s.diff(Pt[i], t, 2).subs(t, 0))
    for i in range(4)
])
print("STAGE 3 COMPLETE", flush=True)

p0 = (gu - gqdot + Pddot).applyfunc(clean)

assert all(x.is_finite for x in p0)

# Exact algebraic decomposition of the computed covector.
p0raised = eta*p0
p0T = clean((p0.T*T)[0])
p0N = clean((p0.T*N)[0])
transverse = (proj*p0raised).applyfunc(clean)

reconstructed = -p0T*T + p0N*N + transverse
residual = (reconstructed-p0raised).applyfunc(clean)
assert residual == s.zeros(4, 1)

print()
print("PARTIAL-u GRADIENT =", gu)
print("q-GRADIENT DERIVATIVE =", gqdot)
print("SECOND P DERIVATIVE =", Pddot)
print("FULL GEOMETRIC p0 =", p0)
print("p0.T =", p0T)
print("p0.N =", p0N)
print("TRANSVERSE p0 =", transverse)
print("RECONSTRUCTION RESIDUAL =", residual)

print()
print("R22 COMPLETE")
print("GEOMETRIC p0: EVALUATED AT ONE FIFTH-ORDER JET")
print("SECOND INDEPENDENT JET: PENDING")
print("HAMILTONIAN RECONSTRUCTION: OPEN")
print("FULL DIRAC SYSTEM: OPEN")
print("DEFI TORSION IDENTIFICATION: UNRESOLVED")
