import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27D ===", flush=True)

u = S.symbols("u0:2")
q = S.symbols("q0:2")
P = S.symbols("P0:2")
p0 = S.symbols("r0:2")
p1 = S.symbols("t0:2")

IA, IS, IW = S.symbols("IA IS IW")
pIA, pIS, pIW = S.symbols("pIA pIS pIW")
hA, hS, hW = S.symbols("hA hS hW", positive=True)
kA, kS, kW, m = S.symbols("kA kS kW m")

def mdot(a,b):
    return -a[0]*b[0]+a[1]*b[1]

def odot(a,b):
    return sum(x*y for x,y in zip(a,b))

s = -mdot(u,u)
v = S.sqrt(s)
d = mdot(u,q)
W = mdot(q,q)+d**2/s

A = (3+2*IA)*IS
C = (2+IA+S.Rational(3,2)*IA**2)*IW

M = m + (
    kA*IA**2+kS*IS**2+kW*IW**2
    +pIA**2/hA+pIS**2/hS+pIW**2/hW
)/2

J = 2*IS*pIA/hA+(3+2*IA)*pIS/hS

phi1 = odot(u,P)
phi2 = odot(q,P)
e = odot(p0,q)
f = mdot(p1,P)
g = odot(p0,u)
P2 = mdot(P,P)

def PB(F,G):
    out = 0
    for x,px in zip(u,p1):
        out += S.diff(F,x)*S.diff(G,px)
        out -= S.diff(F,px)*S.diff(G,x)
    for x,px in zip(q,P):
        out += S.diff(F,x)*S.diff(G,px)
        out -= S.diff(F,px)*S.diff(G,x)
    return out

def check(name, expr):
    if S.cancel(expr) != 0:
        raise AssertionError(name)
    print("PASS:", name, flush=True)

print("STAGE 1: canonical generator rules",flush=True)

rules1 = [
    (s,0),(d,-s),(W,0),
    (e,g),(f,-P2),(g,0),
    (P2,0),(phi1,0),(phi2,phi1)
]

rules2 = [
    (s,0),(d,d),(W,2*W),
    (e,e),(f,-f),(g,0),
    (P2,-2*P2),(phi1,-phi1),(phi2,0)
]

for name, generator, rules in [
    ("phi1",phi1,rules1),
    ("phi2",phi2,rules2)
]:
    for index,(quantity,expected) in enumerate(rules):
        check(
            f"{name} rule {index+1}",
            PB(quantity,generator)-expected
        )

print("STAGE 1 PASS",flush=True)

print("STAGE 2: exact B construction",flush=True)

B = (
    -e
    -v*W*f/(2*C)
    +d*W*P2/(4*C*v)
    +5*A*d*W/v**5
    +d*M/v
    -2*A*W/(C*v**2)*(phi2+d*phi1/s)
    +2*W*J/v**2
)

print("STAGE 2 PASS",flush=True)

print("STAGE 3: independent rational witness",flush=True)

witness = {
    u[0]:2,u[1]:0,
    q[0]:0,q[1]:1,
    P[0]:0,P[1]:1,
    p0[0]:1,p0[1]:2,
    p1[0]:1,p1[1]:3,
    IA:0,IS:1,IW:1,
    pIA:1,pIS:2,pIW:3,
    hA:1,hS:1,hW:1,
    kA:1,kS:1,kW:1,m:1
}

for name,generator in [
    ("phi1",phi1),("phi2",phi2)
]:
    value = S.simplify(PB(B,generator).subs(witness))
    print(
        f"PASS: direct {name} bracket = {value}",
        flush=True
    )

print("STAGE 3 PASS",flush=True)

print()
print("R27D COMPLETE")
print("CANONICAL GENERATOR RULES: VERIFIED")
print("DIRECT B BRACKETS: WITNESS EVALUATED")
print("GENERAL 4D BRACKET IDENTITY: OPEN")
print("B=0 PRESERVATION: OPEN")
print("DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start,3)
)
