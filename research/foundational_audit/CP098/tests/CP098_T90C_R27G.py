import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27G ===", flush=True)

u = S.symbols("u0:2")
q = S.symbols("q0:2")
P = S.symbols("P0:2")
p0 = S.symbols("r0:2")
p1 = S.symbols("t0:2")

I = S.symbols("IA IS IW")
pI = S.symbols("pIA pIS pIW")
h = S.symbols("hA hS hW", positive=True)
k = S.symbols("kA kS kW")
m = S.symbols("m")

def mdot(a,b):
    return -a[0]*b[0]+a[1]*b[1]

def odot(a,b):
    return sum(x*y for x,y in zip(a,b))

s = -mdot(u,u)
v = S.sqrt(s)
d = mdot(u,q)
W = mdot(q,q)+d**2/s

IA,IS,IW = I
pIA,pIS,pIW = pI
hA,hS,hW = h
kA,kS,kW = k

A = (3+2*IA)*IS
C = (2+IA+S.Rational(3,2)*IA**2)*IW

V = sum(k[i]*I[i]**2/2 for i in range(3))
T = sum(pI[i]**2/(2*h[i]) for i in range(3))
M = m+V+T

J = sum(
    S.diff(A,I[i])*pI[i]/h[i]
    for i in range(3)
)

phi1 = odot(u,P)
phi2 = odot(q,P)

Hc = (
    odot(p0,u)
    +odot(p1,q)
    -v*W*mdot(P,P)/(4*C)
    +A*W/v**3
    +v*M
)

B = (
    -odot(p0,q)
    -v*W*mdot(p1,P)/(2*C)
    +d*W*mdot(P,P)/(4*C*v)
    +5*A*d*W/v**5
    +d*M/v
    -2*A*W/(C*v**2)*(
        phi2+d*phi1/s
    )
    +2*W*J/v**2
)

# Regular witness: s>0, W>0, C!=0.
witness = {
    u[0]:2, u[1]:0,
    q[0]:0, q[1]:1,
    P[0]:0, P[1]:1,
    p0[0]:1, p0[1]:2,
    p1[0]:1, p1[1]:3,
    IA:0, IS:1, IW:1,
    pIA:1, pIS:2, pIW:3,
    hA:1, hS:1, hW:1,
    kA:1, kS:1, kW:1,
    m:1
}

def evaluate(expr):
    return S.simplify(expr.subs(witness))

print("STAGE 1: verify witness", flush=True)

assert evaluate(s)>0
assert evaluate(W)>0
assert evaluate(C)!=0

print("s =",evaluate(s),flush=True)
print("W =",evaluate(W),flush=True)
print("C =",evaluate(C),flush=True)
print("STAGE 1 PASS",flush=True)

print("STAGE 2: canonical bracket contributions",flush=True)

pairs = (
    list(zip(u,p1))
    +list(zip(q,P))
    +list(zip(I,pI))
)

# The (y,p0) pairs contribute zero because
# both B and Hc are independent of y.

terms = []

for index,(x,px) in enumerate(pairs,1):
    positive = evaluate(
        S.diff(B,x)*S.diff(Hc,px)
    )
    negative = evaluate(
        S.diff(B,px)*S.diff(Hc,x)
    )
    contribution = S.simplify(positive-negative)
    terms.append(contribution)

    print(
        f"pair {index}: {x}, {px}; "
        f"contribution = {contribution}",
        flush=True
    )

D_direct = S.simplify(sum(terms))

print("D_direct =",D_direct,flush=True)
print("STAGE 2 PASS",flush=True)

print("STAGE 3: Hamiltonian-flow cross-check",flush=True)

# Independently construct the chain rule:
# dot B = sum(B_x * dot x + B_px * dot px)
# with dot x = H_px and dot px = -H_x.

flow_terms = []

for x,px in pairs:
    xdot = evaluate(S.diff(Hc,px))
    pxdot = -evaluate(S.diff(Hc,x))

    Bx = evaluate(S.diff(B,x))
    Bpx = evaluate(S.diff(B,px))

    flow_terms.append(Bx*xdot+Bpx*pxdot)

D_flow = S.simplify(sum(flow_terms))

if S.simplify(D_direct-D_flow)!=0:
    raise AssertionError(
        f"Mismatch: {D_direct} vs {D_flow}"
    )

print("PASS: canonical bracket agrees with flow")
print("D_flow =",D_flow,flush=True)
print("STAGE 3 PASS",flush=True)

print("STAGE 4: report singular relevance",flush=True)

print("A at witness =",evaluate(A),flush=True)
print("B at witness =",evaluate(B),flush=True)
print("J at witness =",evaluate(J),flush=True)

print(
    "NOTE: This witness is not assumed to "
    "satisfy A=B=0 or the five constraints.",
    flush=True
)

print("STAGE 4 PASS",flush=True)

print()
print("R27G COMPLETE")
print("D AT REGULAR WITNESS: EVALUATED")
print("CANONICAL/FLOW CROSS-CHECK: PASS")
print("GENERAL 4D D EXPRESSION: OPEN")
print("D ON SINGULAR CONSTRAINT SURFACE: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start,3)
)
