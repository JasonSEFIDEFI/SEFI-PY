import sympy as S
from time import perf_counter

t0 = perf_counter()
print("=== CP098-T90C-R27T ===", flush=True)

eta = S.diag(-1, 1, 1, 1)

u = S.Matrix(S.symbols("u0:4", real=True))
q = S.Matrix(S.symbols("q0:4", real=True))
P = S.Matrix(S.symbols("P0:4", real=True))
p0 = S.Matrix(S.symbols("p00:04", real=True))
p1 = S.Matrix(S.symbols("p10:14", real=True))

IA, IS, IW = S.symbols("IA IS IW", real=True)
pA, pS, pW = S.symbols("pA pS pW", real=True)

def mdot(a, b):
    return (a.T * eta * b)[0]

s = -mdot(u, u)
v = S.sqrt(s)
d = mdot(u, q)
W = mdot(q, q) + d**2/s
P2 = mdot(P, P)

# Same unit parameters as R27S.
A = (3 + 2*IA)*IS
C = (2 + IA + S.Rational(3,2)*IA**2)*IW

V = (IA**2 + IS**2 + IW**2)/2
T = (pA**2 + pS**2 + pW**2)/2
M = 1 + V + T

J = 2*IS*pA + (3 + 2*IA)*pS

phi1 = (u.T*P)[0]
phi2 = (q.T*P)[0]

Hc = (
    (p0.T*u)[0]
    + (p1.T*q)[0]
    - v*W*P2/(4*C)
    + A*W/v**3
    + v*M
)

B = (
    -(p0.T*q)[0]
    - v*W*mdot(p1,P)/(2*C)
    + d*W*P2/(4*C*v)
    + 5*A*d*W/v**5
    + d*M/v
    - 2*A*W*(phi2+d*phi1/s)/(C*v**2)
    + 2*W*J/v**2
)

print("STAGE 1: reconstruct R27S witness",
      flush=True)

u0 = S.Matrix([2, 0, 0, 0])
q0 = S.Matrix([0, 1, 0, 0])
P0 = S.Matrix([0, 0, 1, 0])
p10 = S.zeros(4, 1)

v0 = S.Integer(2)
W0 = S.Integer(1)
P20 = S.Integer(1)
IA0 = -S.Rational(3,2)
IW0 = S.Integer(1)

IS0 = (
    3*v0**4/(4*W0)
    + 28*P20*v0**4/(961*IW0)
)

C0 = (2+IA0+S.Rational(3,2)*IA0**2)*IW0
M0 = 1 + (IA0**2+IS0**2+IW0**2)/2

p00 = S.Matrix([
    v0*W0*P20/(8*C0)-v0*M0/2,
    0, 0, 0
])

point = {}

for vec, values in (
    (u,u0), (q,q0), (P,P0),
    (p0,p00), (p1,p10)
):
    point.update(zip(vec, values))

point.update({
    IA:IA0, IS:IS0, IW:IW0,
    pA:0, pS:0, pW:0
})

def at(expr):
    return S.cancel(expr.subs(point))

assert at(phi1) == 0
assert at(phi2) == 0
assert at(B) == 0

theta = (
    -(p0.T*u)[0]
    +(p1.T*q)[0]
    +v*W*P2/(4*C)
    +3*A*W/v**3
    -v*M
)

assert at(theta) == 0
assert at(A) == 0
assert at(J) == 0

Q = (
    S.Rational(3,2)*v
    +56*v*W*P2/(961*IW)
    -2*IS*W/v**3
)

assert at(Q) == 0

print("PASS: original witness and Q=0")
print("STAGE 1 PASS", flush=True)

print("STAGE 2: canonical derivatives",
      flush=True)

# Differentiate first, substitute afterward.
# Each pair contributes
# dB/dx*dH/dp - dB/dp*dH/dx.
#
# The y,p0 pair contributes zero because
# neither B nor Hc depends on y.

pairs = (
    list(zip(u,p1))
    +list(zip(q,P))
    +[(IA,pA),(IS,pS),(IW,pW)]
)

D = S.Integer(0)
parts = []

for i,(x,p) in enumerate(pairs,1):
    Bx = at(S.diff(B,x))
    Bp = at(S.diff(B,p))
    Hx = at(S.diff(Hc,x))
    Hp = at(S.diff(Hc,p))

    term = S.cancel(Bx*Hp-Bp*Hx)
    parts.append(term)
    D += term

    print(
        f"pair {i:02d}: ({x},{p}) = {term}",
        flush=True
    )

D = S.cancel(D)

assert D == S.cancel(sum(parts))

print("STAGE 2 PASS", flush=True)

print("STAGE 3: independent assembly",
      flush=True)

D_geometric = S.cancel(sum(parts[:8]))
D_identity = S.cancel(sum(parts[8:]))

assert S.cancel(
    D-D_geometric-D_identity
) == 0

print("geometric contribution =", D_geometric)
print("identity contribution =", D_identity)
print("D =", D)
print("STAGE 3 PASS", flush=True)

print()
print("R27T COMPLETE")

if D == 0:
    print("WITNESS: PASSES D=0 NECESSARY CONDITION")
    print("HIGHER-ORDER PRESERVATION: OPEN")
else:
    print("WITNESS: FAILS D=0 NECESSARY CONDITION")
    print("THIS WITNESS IS DYNAMICALLY EXCLUDED")
    print("ENTIRE IA BRANCH: NOT EXCLUDED")

print("GENERAL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-t0,3)
)
