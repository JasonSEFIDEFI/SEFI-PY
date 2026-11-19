import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R26U ===", flush=True)

def check(name, expression):
    result = S.simplify(expression)
    if result != 0:
        print("FAIL:", name, result, flush=True)
        raise AssertionError(name)
    print("PASS:", name, flush=True)

print("STAGE 1: identity-sector evolution", flush=True)

IA, IS, IW = S.symbols("IA IS IW")
pIA, pIS, pIW = S.symbols("pIA pIS pIW")
hA, hS, hW, v = S.symbols(
    "hA hS hW v", positive=True
)

I = (IA, IS, IW)
pI = (pIA, pIS, pIW)
h = (hA, hS, hW)

A = (3 + 2*IA)*IS

T = sum(pI[i]**2/(2*h[i]) for i in range(3))
H_identity = v*T

Adot_direct = sum(
    S.diff(A,I[i])*S.diff(H_identity,pI[i])
    for i in range(3)
)

Adot_expected = v*sum(
    S.diff(A,I[i])*pI[i]/h[i]
    for i in range(3)
)

check("exact identity evolution",
      Adot_direct-Adot_expected)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: exact constrained witness", flush=True)

# Metric diag(-1,1,1,1).
# u=(1,0,0,0), q=(0,1,0,0).
# Therefore s=v=W=1 and d=0.
#
# mu1=3+2IA
# mu2=2+IA+3IA²/2
#
# IA=IS=0, IW=1.
# h=(1,1,1), k=(0,0,0), m=1.
#
# pI=(0,1,0), p1=P=0.
# p0=(-3/2,6,0,0).

eta = S.diag(-1,1,1,1)

u = S.Matrix([1,0,0,0])
q = S.Matrix([0,1,0,0])
P = S.zeros(4,1)
p1 = S.zeros(4,1)
p0 = S.Matrix([S.Rational(-3,2),6,0,0])

def ordinary(x,y):
    return (x.T*y)[0]

def mdot(x,y):
    return (x.T*eta*y)[0]

s = -mdot(u,u)
speed = S.sqrt(s)
d = mdot(u,q)
W = mdot(q,q)+d**2/s

values = {
    IA:0, IS:0, IW:1,
    pIA:0, pIS:1, pIW:0,
    hA:1, hS:1, hW:1,
    v:1
}

Aw = A.subs(values)
Cw = S.Integer(2)
Mw = S.Rational(3,2)

phi1 = ordinary(u,P)
phi2 = ordinary(q,P)
chi = ordinary(u,p1)

psi = ordinary(q,p1)+2*Aw*W/speed**3

Hc = (
    ordinary(p0,u)
    + ordinary(p1,q)
    - speed*W*mdot(P,P)/(4*Cw)
    + Aw*W/speed**3
    + speed*Mw
)

theta = 2*psi-Hc

check("s=1",s-1)
check("W=1",W-1)
check("A=0",Aw)
check("C=2",Cw-2)

for name, expression in (
    ("phi1",phi1),
    ("phi2",phi2),
    ("chi",chi),
    ("psi",psi),
    ("theta",theta)
):
    check(name+"=0",expression)

print("STAGE 2 PASS", flush=True)

print("STAGE 3: B and tangency", flush=True)

# At this witness:
# d=0, A=0, p1=P=0.
#
# General B reduces to:
# -p0.q + 2W/v² sum(A_i*pI_i/h_i).

identity_sum = sum(
    (S.diff(A,I[i])*pI[i]/h[i]).subs(values)
    for i in range(3)
)

Bw = (
    -ordinary(p0,q)
    +2*W/speed**2*identity_sum
)

Adotw = Adot_expected.subs(values)

check("B=0",Bw)
check("Adot=3",Adotw-3)

print("WITNESS B =",Bw)
print("WITNESS Adot =",Adotw)

print("STAGE 3 PASS", flush=True)

print()
print("R26U COMPLETE")
print("EXACT FIVE-CONSTRAINT WITNESS: PASS")
print("A=0 AND B=0: VERIFIED")
print("Adot NONZERO: VERIFIED")
print("A=0,B=0 NOT AUTOMATICALLY INVARIANT")
print("INVARIANT SINGULAR SUBBRANCH: OPEN")
print("ADDITIONAL CONSISTENCY CONDITIONS: OPEN")
print("GAUGE GENERATORS: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(monotonic()-start,3))
