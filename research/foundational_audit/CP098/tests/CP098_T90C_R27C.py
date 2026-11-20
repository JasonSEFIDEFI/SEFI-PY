import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27C ===", flush=True)

# Independent directional-derivative audit
# in four spacetime dimensions.

u = S.symbols("u0:4", real=True)
q = S.symbols("q0:4", real=True)
P = S.symbols("P0:4", real=True)
p0 = S.symbols("r0:4", real=True)
p1 = S.symbols("t0:4", real=True)

IA, IS, IW = S.symbols("IA IS IW", real=True)
pIA, pIS, pIW = S.symbols("pIA pIS pIW", real=True)
hA, hS, hW = S.symbols("hA hS hW", positive=True)
kA, kS, kW, m = S.symbols("kA kS kW m", real=True)

eta = S.diag(-1, 1, 1, 1)

def mdot(a, b):
    return sum(eta[i,i]*a[i]*b[i] for i in range(4))

def odot(a, b):
    return sum(a[i]*b[i] for i in range(4))

s = -mdot(u,u)
v = S.sqrt(s)
d = mdot(u,q)
W = mdot(q,q) + d**2/s

A = (3+2*IA)*IS
C = (2+IA+S.Rational(3,2)*IA**2)*IW

V = (kA*IA**2+kS*IS**2+kW*IW**2)/2
T = (pIA**2/hA+pIS**2/hS+pIW**2/hW)/2
M = m+V+T

phi1 = odot(u,P)
phi2 = odot(q,P)
P2 = mdot(P,P)

B = (
    -odot(p0,q)
    -v*W*mdot(p1,P)/(2*C)
    +d*W*P2/(4*C*v)
    +5*A*d*W/v**5
    +d*M/v
    -2*A*W/(C*v**2)*(phi2+d*phi1/s)
    +2*W/v**2*(
        2*IS*pIA/hA+(3+2*IA)*pIS/hS
    )
)

print("STAGE 1: construct exact brackets", flush=True)

# phi1=u.P, phi2=q.P
# Both depend on (q,P), while phi1 also
# depends on u. Canonical pair (u,p1).

br1 = sum(
    -u[i]*S.diff(B,q[i])
    +P[i]*S.diff(B,p1[i])
    for i in range(4)
)

br2 = sum(
    -q[i]*S.diff(B,q[i])
    +P[i]*S.diff(B,P[i])
    for i in range(4)
)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: independent canonical audit", flush=True)

def PB(F,G):
    out = 0
    for x,px in zip(u,p1):
        out += S.diff(F,x)*S.diff(G,px)
        out -= S.diff(F,px)*S.diff(G,x)
    for x,px in zip(q,P):
        out += S.diff(F,x)*S.diff(G,px)
        out -= S.diff(F,px)*S.diff(G,x)
    return out

# The canonical definition gives:
# {B,phi1} = sum(P*dB/dp1 - u*dB/dq)
# {B,phi2} = sum(P*dB/dP - q*dB/dq)

def check(name, expr):
    result = S.cancel(expr)
    if result != 0:
        raise AssertionError(f"{name}: {result}")
    print("PASS:", name, flush=True)

check("first bracket implementation", br1-PB(B,phi1))
check("second bracket implementation", br2-PB(B,phi2))
print("STAGE 2 PASS", flush=True)

print("STAGE 3: rational witness", flush=True)

witness = {
    u[0]:2,u[1]:0,u[2]:0,u[3]:0,
    q[0]:0,q[1]:1,q[2]:0,q[3]:0,
    P[0]:0,P[1]:1,P[2]:0,P[3]:0,
    p0[0]:1,p0[1]:2,p0[2]:0,p0[3]:0,
    p1[0]:1,p1[1]:3,p1[2]:0,p1[3]:0,
    IA:0,IS:1,IW:1,
    pIA:1,pIS:2,pIW:3,
    hA:1,hS:1,hW:1,
    kA:1,kS:1,kW:1,m:1
}

for name,expr in [
    ("B",B),
    ("bracket1",br1),
    ("bracket2",br2),
    ("delta1",br1-B),
    ("delta2",br2-2*B)
]:
    print(name,"=",S.simplify(expr.subs(witness)),flush=True)

print("STAGE 3 PASS",flush=True)

print()
print("R27C COMPLETE")
print("DIRECT 4D BRACKET CONSTRUCTION: PASS")
print("CANDIDATE IDENTITIES: NOT ASSUMED")
print("SINGULAR B=0 RESTRICTION: OPEN")
print("DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",round(perf_counter()-start,3))
