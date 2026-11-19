import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R26O ===", flush=True)

eta = S.diag(-1, 1, 1, 1)

def pair(a, b):
    return (a.T*b)[0]

def mdot(a, b):
    return (a.T*eta*b)[0]

u = S.Matrix(S.symbols("u0:4"))
q = S.Matrix(S.symbols("q0:4"))
P = S.Matrix(S.symbols("P0:4"))
p0 = S.Matrix(S.symbols("p0_0:4"))
p1 = S.Matrix(S.symbols("p1_0:4"))

IA, IS, IW = S.symbols("IA IS IW")
pa, ps, pw = S.symbols("pa ps pw")
m = S.symbols("m")

I = [IA, IS, IW]
pI = [pa, ps, pw]

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

Hc = (
    pair(p0,u)+pair(p1,q)
    -v*W*P2/(4*C)
    +A*W/v**3+v*M
)

theta = (
    -pair(p0,u)+pair(p1,q)
    +v*W*P2/(4*C)
    +3*A*W/v**3-v*M
)

names = ["phi1","phi2","chi","psi","theta"]
F = [phi1,phi2,chi,psi,theta]

canonical = (
    list(zip(u,p1))
    +list(zip(q,P))
    +list(zip(I,pI))
)

print("STAGE 1: exact definitions",flush=True)

assert S.expand(Hc+theta-2*psi)==0
assert len(canonical)==11
assert len(F)==5

print("Hc+theta=2psi: PASS")
print("CANONICAL PAIRS =",len(canonical))
print("STAGE 1 PASS",flush=True)

print("STAGE 2: prepare independent derivatives",
      flush=True)

# Compute gradients of the actual constraint
# functions, without using the bracket inventory.
#
# Do not expand the complete Poisson matrix
# symbolically. Evaluate gradients at each jet.

gradients = []

for i,f in enumerate(F):
    dx = []
    dp = []

    for x,p in canonical:
        dx.append(S.diff(f,x))
        dp.append(S.diff(f,p))

    gradients.append((dx,dp))
    print("GRADIENT",names[i],"READY",flush=True)

print("STAGE 2 PASS",flush=True)

print("STAGE 3: construct predicted bracket",
      flush=True)

identity_flow = sum(
    S.diff(A,x)*p
    for x,p in zip(I,pI)
)

B = (
    -pair(p0,q)
    -v*W*pair(p1,eta*P)/(2*C)
    +d*W*P2/(4*C*v)
    +5*A*d*W/v**5
    +d*M/v
    -2*A*W*(phi2+d*phi1/s)/(C*v**2)
    +2*W*identity_flow/v**2
)

b = 2*A*W/v**3

# This matrix is reconstructed from the
# earlier bracket inventory, independently
# of the derivatives calculated above.

expected = S.zeros(5)

entries = {
    (0,1):-phi1,
    (0,2):phi1,
    (0,3):phi2-chi,
    (0,4):phi2-chi,
    (1,2):S.Integer(0),
    (1,3):-psi-b,
    (1,4):-psi-2*b,
    (2,3):psi+2*b,
    (2,4):2*psi+4*b-theta,
    (3,4):-B
}

assert len(entries)==10

for (i,j),value in entries.items():
    expected[i,j] = value
    expected[j,i] = -value

assert expected+expected.T==S.zeros(5)

print("TEN EXPECTED BRACKETS READY")
print("STAGE 3 PASS",flush=True)

jets = [
    {
        "name":"FIVE-CONSTRAINT WITNESS",
        "u":[5,3,0,0],
        "q":[1,2,2,1],
        "P":[0,0,0,0],
        "p0":[-S.Rational(52861,5120),0,0,0],
        "p1":[0,0,-S.Rational(387,1024),0],
        "I":[0,1,1],
        "pI":[2,-1,3],
        "m":5,
        "constrained":True
    },
    {
        "name":"OFF-SURFACE WITNESS",
        "u":[13,3,4,0],
        "q":[2,3,1,4],
        "P":[1,2,-1,3],
        "p0":[2,-1,3,1],
        "p1":[1,0,2,-2],
        "I":[1,2,1],
        "pI":[1,3,-2],
        "m":7,
        "constrained":False
    }
]

def exact_at(expr,point):
    return S.simplify(expr.subs(point))

def check_jet(jet):
    print("CHECK:",jet["name"],flush=True)

    point = {}

    for variables,key in [
        (u,"u"),(q,"q"),(P,"P"),
        (p0,"p0"),(p1,"p1"),
        (I,"I"),(pI,"pI")
    ]:
        point.update(zip(variables,jet[key]))

    point[m] = jet["m"]

    assert exact_at(s,point)>0
    assert exact_at(W,point)>0
    assert exact_at(C,point)!=0
    assert exact_at(A,point)!=0

    values = [
        exact_at(f,point)
        for f in F
    ]

    if jet["constrained"]:
        assert all(x==0 for x in values)
        print("FIVE CONSTRAINTS: ZERO",flush=True)
    else:
        assert any(x!=0 for x in values)
        print("OFF-SURFACE: CONFIRMED",flush=True)

    evaluated = []

    for i,(dx,dp) in enumerate(gradients):
        gx = [exact_at(z,point) for z in dx]
        gp = [exact_at(z,point) for z in dp]
        evaluated.append((gx,gp))

        print(
            "EVALUATED GRADIENT",
            names[i],
            flush=True
        )

    actual = S.zeros(5)

    for i in range(5):
        for j in range(i+1,5):
            xi,pi = evaluated[i]
            xj,pj = evaluated[j]

            value = sum(
                xi[k]*pj[k]-pi[k]*xj[k]
                for k in range(11)
            )

            actual[i,j] = S.factor(value)
            actual[j,i] = -actual[i,j]

    target = expected.subs(point)

    for i in range(5):
        for j in range(i+1,5):
            residual = S.factor(
                actual[i,j]-target[i,j]
            )

            print(
                names[i],names[j],
                "RESIDUAL =",residual,
                flush=True
            )

            assert residual==0, (
                jet["name"],
                names[i],names[j],
                residual
            )

    assert actual+actual.T==S.zeros(5)

    print("ALL TEN BRACKETS: PASS",flush=True)

    if jet["constrained"]:
        rank = actual.rank()
        assert rank==2
        assert len(actual.nullspace())==3

        print("WEAK RANK =",rank)
        print("WEAK NULLITY = 3")

    return point

print("STAGE 4: constrained canonical audit",
      flush=True)

check_jet(jets[0])
print("STAGE 4 PASS",flush=True)

print("STAGE 5: off-surface canonical audit",
      flush=True)

check_jet(jets[1])
print("STAGE 5 PASS",flush=True)

print()
print("R26O COMPLETE")
print("TEN CANONICAL BRACKETS: DIRECTLY CHECKED")
print("FIVE-CONSTRAINT WITNESS: PASS")
print("OFF-SURFACE WITNESS: PASS")
print("WEAK MATRIX RANK: 2 AT WITNESS")
print("WEAK NULLITY: 3 AT WITNESS")
print("GENERAL SYMBOLIC MATRIX: NOT PROVEN HERE")
print("A-NONZERO CLASSIFICATION: STRONGLY SUPPORTED")
print("A-ZERO BRANCH: OPEN")
print("GAUGE GENERATORS: OPEN")
print("DEFI CORRESPONDENCE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(monotonic()-start,3)
)
