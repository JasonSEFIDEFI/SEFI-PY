import sympy as S
from time import monotonic

t0 = monotonic()
print("=== CP098-T90C-R26I ===", flush=True)

eta = S.diag(-1, 1, 1, 1)

def pair(a,b):
    return (a.T*b)[0]

def mdot(a,b):
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

point = {
    u[0]:5,u[1]:3,u[2]:0,u[3]:0,
    q[0]:1,q[1]:2,q[2]:2,q[3]:1,
    IA:0,IS:1,IW:1,
    pa:2,ps:-1,pw:3,
    m:5
}

for x in list(P)+list(p0)+list(p1):
    point[x] = S.Integer(0)

point[p1[2]] = -S.Rational(387,1024)

def at(expr):
    return S.simplify(expr.subs(point))

print("STAGE 1: five-constraint witness", flush=True)

# Choose p0 so that theta=0 exactly.
# At this point u=(5,3,0,0).
# Hence p0.u=5*p0[0] when other
# p0 components are zero.

base = at(
    pair(p1,q)
    +v*W*P2/(4*C)
    +3*A*W/v**3
    -v*M
)

point[p0[0]] = base/5

assert at(v)==4
assert at(W)==S.Rational(129,16)
assert at(A)==3
assert at(C)==2

for name, expr in [
    ("phi1",phi1),
    ("phi2",phi2),
    ("chi",chi),
    ("psi",psi),
    ("theta",theta)
]:
    value = at(expr)
    print(name, "=", value, flush=True)
    assert value == 0

print("p0[0] =",point[p0[0]])
print("STAGE 1 PASS", flush=True)

print("STAGE 2: exact canonical derivatives",
      flush=True)

# H and theta have no y dependence.
# The (y,p0) bracket therefore vanishes.
# Remaining canonical pairs are:
# (u,p1), (q,P), and (I,pI).

canonical = (
    list(zip(u,p1))
    +list(zip(q,P))
    +list(zip(I,pI))
)

terms = []

for index,(x,p) in enumerate(canonical,1):
    term = (
        S.diff(theta,x)*S.diff(H,p)
        -S.diff(theta,p)*S.diff(H,x)
    )
    value = at(term)
    terms.append(value)
    print(
        "PAIR",index,"/",len(canonical),
        "DONE",
        flush=True
    )

bracket = S.factor(sum(terms))

print("{theta,Hc} =",bracket,flush=True)
print("STAGE 2 PASS",flush=True)

print("STAGE 3: Jacobi multiplier brackets",
      flush=True)

# Previously established exact identities:
#
# theta={chi,Hc}
# {chi,phi1}=-phi1
# {chi,phi2}=0
# {Hc,phi1}=chi-phi2
# {Hc,phi2}=psi
# {chi,psi}=psi+4AW/v^3
#
# Jacobi implies:
#
# {theta,phi1}=chi-phi2
# {theta,phi2}=psi+4AW/v^3.

b1 = at(chi-phi2)
b2 = at(psi+4*A*W/v**3)

print("{theta,phi1} WEAK =",b1)
print("{theta,phi2} WEAK =",b2)

assert b1 == 0
assert b2 == S.Rational(387,256)

print("STAGE 3 PASS",flush=True)

print("STAGE 4: psi-fixed multiplier",flush=True)

# From R26H at the corresponding
# four-constraint witness:
lambda2 = -S.Rational(71695,12384)

# theta preservation:
# dot(theta)={theta,Hc}
#            +lambda1*{theta,phi1}
#            +lambda2*{theta,phi2}.

remainder = S.factor(
    bracket + lambda2*b2
)

print("lambda2 =",lambda2)
print("dot(theta) WEAK =",remainder)

if remainder == 0:
    print(
        "WITNESS: theta preservation satisfied",
        flush=True
    )
else:
    print(
        "WITNESS: theta preservation imposes",
        "an additional condition",
        flush=True
    )

print("STAGE 4 PASS",flush=True)

print()
print("R26I COMPLETE")
print("FIVE-CONSTRAINT WITNESS: VERIFIED")
print("DIRECT theta HAMILTONIAN BRACKET: EVALUATED")
print("MULTIPLIER DEPENDENCE: EVALUATED")
print("GENERAL TERTIARY PRESERVATION: OPEN")
print("ADDITIONAL CONSTRAINT CLASSIFICATION: OPEN")
print("FULL DIRAC CLASSIFICATION: OPEN")
print("PHYSICAL DEGREES OF FREEDOM: NOT COUNTED")
print("DEFI CORRESPONDENCE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",round(monotonic()-t0,3))
