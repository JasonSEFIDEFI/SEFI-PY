import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R26D ===", flush=True)

eta = S.diag(-1, 1, 1, 1)

def dot(a, b):
    return (a.T * eta * b)[0]

def pair(a, b):
    return (a.T * b)[0]

u = S.Matrix(S.symbols("u0:4", real=True))
q = S.Matrix(S.symbols("q0:4", real=True))
P = S.Matrix(S.symbols("P0:4", real=True))

p0 = S.Matrix(S.symbols("p0_0:4", real=True))
p1 = S.Matrix(S.symbols("p1_0:4", real=True))

A, C, m = S.symbols("A C m", real=True)
V, T = S.symbols("V T", real=True)

# A,C,V,T are independent of u,q,P,p0,p1
# for the brackets calculated in this test.
#
# T = sum_a pIa**2/(2*ha).
# The identity variables remain fixed under
# phi1 and phi2.

s = -dot(u, u)
v = S.sqrt(s)
d = dot(u, q)
W = S.expand(dot(q, q) + d**2/s)

P2 = dot(P, P)

phi1 = pair(u, P)
phi2 = pair(q, P)

print("STAGE 1: local Hamiltonian", flush=True)

Hc = (
    pair(p0, u)
    + pair(p1, q)
    - v*W*P2/(4*C)
    + A*W/v**3
    + m*v
    + v*(V + T)
)

print("LOCAL DOMAIN: s>0, W>0, C!=0")
print("Hc CONSTRUCTED")
print("STAGE 1 PASS", flush=True)

# Use only canonical pairs relevant to these
# primary-constraint brackets:
# (u,p1) and (q,P).
#
# The omitted pairs contribute zero because
# phi1 and phi2 do not depend on them.

def PB(f, g):
    total = S.Integer(0)
    for x, p in (
        list(zip(u, p1)) + list(zip(q, P))
    ):
        total += (
            S.diff(f, x)*S.diff(g, p)
            - S.diff(f, p)*S.diff(g, x)
        )
    return S.expand(total)

print("STAGE 2: primary bracket cross-check",
      flush=True)

b12 = PB(phi1, phi2)
b21 = PB(phi2, phi1)

assert S.simplify(b12 + phi1) == 0
assert S.simplify(b21 - phi1) == 0

print("{phi1,phi2} + phi1 =",
      S.simplify(b12 + phi1))
print("{phi2,phi1} - phi1 =",
      S.simplify(b21 - phi1))

print("STAGE 2 PASS", flush=True)

print("STAGE 3: phi1 preservation", flush=True)

pres1 = PB(phi1, Hc)
res1 = S.simplify(pres1 + phi2)

print("{phi1,Hc} + phi2 =", res1)

assert res1 == 0

print("{phi1,Hc} = -phi2")
print("STAGE 3 PASS", flush=True)

print("STAGE 4: phi2 preservation", flush=True)

pres2 = PB(phi2, Hc)

psi = pair(p1, q) + 2*A*W/v**3

res2 = S.simplify(pres2 - psi)

print("{phi2,Hc} - psi =", res2)

assert res2 == 0

print("SECONDARY CANDIDATE:")
print("psi = p1.q + 2*A*W/v**3")

print("STAGE 4 PASS", flush=True)

print("STAGE 5: total-Hamiltonian equations",
      flush=True)

lambda1, lambda2 = S.symbols(
    "lambda1 lambda2", real=True
)

HT = Hc + lambda1*phi1 + lambda2*phi2

dotphi1 = PB(phi1, HT)
dotphi2 = PB(phi2, HT)

expected1 = -phi2 - lambda2*phi1
expected2 = psi + lambda1*phi1

check1 = S.simplify(dotphi1 - expected1)
check2 = S.simplify(dotphi2 - expected2)

print("dot(phi1) RESIDUAL =", check1)
print("dot(phi2) RESIDUAL =", check2)

assert check1 == 0
assert check2 == 0

print("ON PRIMARY SURFACE:")
print("dot(phi1) = 0")
print("dot(phi2) = psi")
print("STAGE 5 PASS", flush=True)

print()
print("R26D COMPLETE")
print("LOCAL CANONICAL HAMILTONIAN: CONSTRUCTED")
print("PRIMARY BRACKET CROSS-CHECK: PASS")
print("phi1 PRESERVATION: PASS")
print("phi2 PRESERVATION: PASS")
print("SECONDARY CONSTRAINT CANDIDATE: DERIVED")
print("SECONDARY PRESERVATION: OPEN")
print("FULL CONSTRAINT CLASSIFICATION: OPEN")
print("DIRAC ALGORITHM: INCOMPLETE")
print("DEFI CORRESPONDENCE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(monotonic()-start, 3))
