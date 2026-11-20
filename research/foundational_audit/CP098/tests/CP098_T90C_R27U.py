import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27U ===", flush=True)

t = S.symbols("t", real=True)
eta = S.diag(-1, 1, 1, 1)

u = S.Matrix(S.symbols("u0:4", real=True))
q = S.Matrix(S.symbols("q0:4", real=True))
P = S.Matrix(S.symbols("P0:4", real=True))
p0 = S.Matrix(S.symbols("p00:04", real=True))
p1 = S.Matrix(S.symbols("p10:14", real=True))

IA, IS, IW = S.symbols("IA IS IW", real=True)
pA, pS, pW = S.symbols("pA pS pW", real=True)

def dot(a, b):
    return (a.T * eta * b)[0]

s = -dot(u, u)
v = S.sqrt(s)
d = dot(u, q)
W = dot(q, q) + d**2/s
P2 = dot(P, P)

A = (3+2*IA)*IS
C = (2+IA+S.Rational(3,2)*IA**2)*IW
M = 1+(IA**2+IS**2+IW**2
       +pA**2+pS**2+pW**2)/2
J = 2*IS*pA+(3+2*IA)*pS

phi1 = (u.T*P)[0]
phi2 = (q.T*P)[0]

H = (
    (p0.T*u)[0]+(p1.T*q)[0]
    -v*W*P2/(4*C)
    +A*W/v**3+v*M
)

B = (
    -(p0.T*q)[0]
    -v*W*dot(p1,P)/(2*C)
    +d*W*P2/(4*C*v)
    +5*A*d*W/v**5
    +d*M/v
    -2*A*W*(phi2+d*phi1/s)/(C*v**2)
    +2*W*J/v**2
)

chi = (u.T*p1)[0]
psi = (q.T*p1)[0]+2*A*W/v**3
theta = (
    -(p0.T*u)[0]+(p1.T*q)[0]
    +v*W*P2/(4*C)
    +3*A*W/v**3-v*M
)

Q = (
    S.Rational(3,2)*v
    +56*v*W*P2/(961*IW)
    -2*IS*W/v**3
)

u0 = S.Matrix([2,0,0,0])
q0 = S.Matrix([0,1,0,0])
P0 = S.Matrix([0,0,1,0])
p10 = S.Matrix([0,0,0,t])

IS0 = S.Rational(11980,961)
C0 = S.Rational(31,8)
M0 = 1+(S.Rational(9,4)+IS0**2+1)/2

p00 = S.Matrix([
    S.Rational(1,4)/C0-M0,
    0,0,0
])

point = {}
for vec, vals in (
    (u,u0),(q,q0),(P,P0),
    (p0,p00),(p1,p10)
):
    point.update(zip(vec,vals))

point.update({
    IA:-S.Rational(3,2),
    IS:IS0,IW:1,pA:0,pS:0,pW:0
})

def at(expr):
    return S.cancel(expr.subs(point))

print("STAGE 1: family constraint audit",
      flush=True)

for name,expr in [
    ("phi1",phi1),("phi2",phi2),
    ("chi",chi),("psi",psi),
    ("theta",theta),("B",B),
    ("A",A),("J",J),("Q",Q),
    ("pA",pA)
]:
    result = at(expr)
    assert result == 0, (name,result)
    print("PASS:",name)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: exact D(t)", flush=True)

pairs = (
    list(zip(u,p1))
    +list(zip(q,P))
    +[(IA,pA),(IS,pS),(IW,pW)]
)

parts = []
for x,p in pairs:
    term = S.cancel(
        at(S.diff(B,x))*at(S.diff(H,p))
        -at(S.diff(B,p))*at(S.diff(H,x))
    )
    parts.append(term)

D = S.factor(sum(parts))

baseline = S.Rational(593952197,14776336)
assert S.cancel(D.subs(t,0)-baseline) == 0

print("PASS: R27T baseline reproduced")
print("D(t) =",D)
print("degree =",S.Poly(D,t).degree())
print("STAGE 2 PASS", flush=True)

print("STAGE 3: real-root analysis",
      flush=True)

roots = S.solveset(
    S.Eq(D,0),t,domain=S.S.Reals
)

print("real solutions =",roots)

if roots == S.S.EmptySet:
    print("NO REAL D=0 POINT IN THIS FAMILY")
else:
    print("D=0 CANDIDATES: REQUIRE FURTHER AUDIT")

print("STAGE 3 PASS", flush=True)

print()
print("R27U COMPLETE")
print("ONE-PARAMETER CONSTRAINT FAMILY: VERIFIED")
print("D=0 ROOTS: REPORTED")
print("GENERAL BRANCH CONSISTENCY: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(perf_counter()-start,3))
