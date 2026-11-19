import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R26J2 ===", flush=True)

eta = S.diag(-1,1,1,1)

def pair(a,b):
    return (a.T*b)[0]

def mdot(a,b):
    return (a.T*eta*b)[0]

u = S.Matrix(S.symbols("u0:4"))
q = S.Matrix(S.symbols("q0:4"))
P = S.Matrix(S.symbols("P0:4"))
p0 = S.Matrix(S.symbols("p0_0:4"))
p1 = S.Matrix(S.symbols("p1_0:4"))

IA,IS,IW = S.symbols("IA IS IW")
pa,ps,pw = S.symbols("pa ps pw")
m = S.symbols("m")

I = [IA,IS,IW]
pI = [pa,ps,pw]

A = (3+2*IA)*IS
C = (2+IA+S.Rational(3,2)*IA**2)*IW

s = -mdot(u,u)
v = S.sqrt(s)
d = mdot(u,q)
W = mdot(q,q)+d**2/s
P2 = mdot(P,P)

V = sum(x*x/2 for x in I)
T = sum(x*x/2 for x in pI)
M = m+V+T

phi1 = pair(u,P)
phi2 = pair(q,P)
chi = pair(u,p1)
psi = pair(q,p1)+2*A*W/v**3

H = (
    pair(p0,u)+pair(p1,q)
    -v*W*P2/(4*C)
    +A*W/v**3+v*M
)

theta = (
    -pair(p0,u)+pair(p1,q)
    +v*W*P2/(4*C)
    +3*A*W/v**3-v*M
)

canonical = (
    list(zip(u,p1))
    +list(zip(q,P))
    +list(zip(I,pI))
)

base_point = {
    u[0]:5,u[1]:3,u[2]:0,u[3]:0,
    q[0]:1,q[1]:2,q[2]:2,q[3]:1,
    IA:0,IS:1,IW:1,
    pa:2,ps:-1,pw:3,
    m:5
}

for x in list(P)+list(p0)+list(p1):
    base_point[x] = S.Integer(0)

base_point[p1[2]] = -S.Rational(387,1024)
base_point[p0[0]] = -S.Rational(52861,5120)

def evaluate(expr,point):
    return S.simplify(expr.subs(point))

def bracket_at(f,g,point):
    total = S.Integer(0)
    for index,(x,p) in enumerate(canonical,1):
        term = (
            S.diff(f,x)*S.diff(g,p)
            -S.diff(f,p)*S.diff(g,x)
        )
        total += evaluate(term,point)
    return S.factor(total)

def run_witness(t):
    point = dict(base_point)

    # Keep p0.u = 5*p0[0]+3*p0[1]
    # constant while changing transverse p0.
    point[p0[1]] = S.Integer(t)
    point[p0[0]] = (
        base_point[p0[0]]
        -S.Rational(3,5)*t
    )

    for name,expr in [
        ("phi1",phi1),
        ("phi2",phi2),
        ("chi",chi),
        ("psi",psi),
        ("theta",theta)
    ]:
        value = evaluate(expr,point)
        assert value == 0, (t,name,value)

    B = bracket_at(psi,H,point)
    D = bracket_at(theta,H,point)

    coefficient = evaluate(2*A*W/v**3,point)
    assert coefficient != 0

    lambda2 = S.factor(-B/coefficient)

    # On the five-constraint surface,
    # {theta,phi2}=4AW/v^3.
    theta_phi2 = evaluate(4*A*W/v**3,point)

    Omega = S.factor(
        D + lambda2*theta_phi2
    )

    print(
        "t =",t,
        "B =",B,
        "D =",D,
        "lambda2 =",lambda2,
        "Omega =",Omega,
        flush=True
    )

    return Omega,B,D

print("STAGE 1: original witness",flush=True)

omega0,B0,D0 = run_witness(0)

assert B0 == S.Rational(1204251,81920)
assert D0 == S.Rational(1204251,40960)
assert omega0 == 0

print("STAGE 1 PASS",flush=True)

print("STAGE 2: transverse momentum variation",
      flush=True)

omega1,B1,D1 = run_witness(1)

delta = S.factor(omega1-omega0)

print("DELTA OMEGA =",delta,flush=True)

if delta == 0:
    print(
        "INDEPENDENCE NOT ESTABLISHED BY THIS",
        "VARIATION",
        flush=True
    )
else:
    print(
        "SIXTH CONDITION VARIES WHILE",
        "FIVE CONSTRAINTS REMAIN ZERO",
        flush=True
    )

print("STAGE 2 PASS",flush=True)

print()
print("R26J2 COMPLETE")
print("TWO FIVE-CONSTRAINT WITNESSES: VERIFIED")
print("lambda2 RECALCULATED AT EACH WITNESS")
print("SIXTH CONDITION VARIATION =",delta)

if delta != 0:
    print("LOCAL INDEPENDENCE EVIDENCE: PASS")
else:
    print("LOCAL INDEPENDENCE EVIDENCE: OPEN")

print("GENERAL SIXTH CONSTRAINT: OPEN")
print("SIXTH PRESERVATION: OPEN")
print("FULL DIRAC CLASSIFICATION: OPEN")
print("DEFI CORRESPONDENCE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(monotonic()-start,3))
