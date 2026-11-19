import sympy as S
from time import monotonic

t0 = monotonic()
print("=== CP098-T90C-R26H ===", flush=True)

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

IA,IS,IW = S.symbols("IA IS IW")
pa,ps,pw = S.symbols("pa ps pw")
m = S.symbols("m")

I = [IA,IS,IW]
pI = [pa,ps,pw]
h = [S.Integer(1)]*3

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

print("STAGE 1: regular constrained witness",
      flush=True)

point = {
    u[0]:5,u[1]:3,u[2]:0,u[3]:0,
    q[0]:1,q[1]:2,q[2]:2,q[3]:1,
    IA:0,IS:1,IW:1,
    pa:2,ps:-1,pw:3,
    m:5
}

for x in P:
    point[x] = S.Integer(0)

for x in p0:
    point[x] = S.Integer(0)

for x in p1:
    point[x] = S.Integer(0)

point[p1[2]] = -S.Rational(387,1024)

def at(expr):
    return S.simplify(expr.subs(point))

assert at(v)==4
assert at(W)==S.Rational(129,16)
assert at(A)==3
assert at(C)==2

for name,expr in [
    ("phi1",phi1),
    ("phi2",phi2),
    ("chi",chi),
    ("psi",psi)
]:
    value=at(expr)
    print(name,"=",value,flush=True)
    assert value==0

print("STAGE 1 PASS",flush=True)

print("STAGE 2: direct canonical bracket",
      flush=True)

# Evaluate each differentiated term at the
# rational witness before summation.
# This prevents full expression expansion.

canonical = (
    list(zip(u,p1))
    +list(zip(q,P))
    +list(zip(I,pI))
)

terms = []

for x,p in canonical:
    term = (
        S.diff(psi,x)*S.diff(H,p)
        -S.diff(psi,p)*S.diff(H,x)
    )
    terms.append(at(term))

direct = S.factor(sum(terms))

print("DIRECT {psi,H} =",direct,flush=True)

print("STAGE 2 PASS",flush=True)

print("STAGE 3: independently derived formula",
      flush=True)

identity_term = sum(
    S.diff(A,x)*p/hh
    for x,p,hh in zip(I,pI,h)
)

predicted = (
    -pair(p0,q)
    +d*W*P2/(4*C*v)
    +5*A*d*W/v**5
    +d*M/v
    -v*W*pair(p1,eta*P)/(2*C)
    -2*A*W*(phi2+d*phi1/s)/(C*v**2)
    +2*W*identity_term/v**2
)

expected = at(predicted)

print("FORMULA VALUE =",expected,flush=True)
print("RESIDUAL =",S.factor(direct-expected))

assert direct==expected

print("STAGE 3 PASS",flush=True)

print("STAGE 4: identity-sector contribution",
      flush=True)

id_value = at(2*W*identity_term/v**2)

print("IDENTITY CONTRIBUTION =",id_value)
assert id_value != 0

without_identity = at(
    predicted-2*W*identity_term/v**2
)

assert direct-without_identity==id_value

print("IDENTITY COUPLING INCLUDED: PASS")
print("STAGE 4 PASS",flush=True)

print("STAGE 5: multiplier equation",
      flush=True)

# Previously verified:
# {psi,phi1}=chi-phi2
# {psi,phi2}=psi+2AW/v^3.
#
# On phi1=phi2=chi=psi=0:
#
# dot(psi)=B+lambda2*(2AW/v^3).
#
# For A!=0, this determines lambda2.

coefficient=at(2*A*W/v**3)
assert coefficient==S.Rational(387,512)
assert coefficient!=0

lambda2_value=S.factor(-direct/coefficient)

print("WEAK MULTIPLIER COEFFICIENT =",
      coefficient)
print("lambda2 AT WITNESS =",lambda2_value)

assert S.factor(
    direct+lambda2_value*coefficient
)==0

print("STAGE 5 PASS",flush=True)

print()
print("R26H COMPLETE")
print("DIRECT CANONICAL BRACKET: EVALUATED")
print("PROPOSED FORMULA: VERIFIED AT WITNESS")
print("IDENTITY CONTRIBUTION: NONZERO")
print("psi PRESERVATION: MULTIPLIER FIXED")
print("GENERAL SYMBOLIC BRACKET: NOT PROVEN")
print("theta PRESERVATION: OPEN")
print("FULL DIRAC CLASSIFICATION: OPEN")
print("DEFI CORRESPONDENCE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(monotonic()-t0,3))
