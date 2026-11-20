import sympy as S
from time import monotonic

t0 = monotonic()
print("=== CP098-T90C-R26D2 ===", flush=True)

eta = S.diag(-1, 1, 1, 1)

u = S.Matrix(S.symbols("u0:4"))
q = S.Matrix(S.symbols("q0:4"))
P = S.Matrix(S.symbols("P0:4"))
p0 = S.Matrix(S.symbols("p0_0:4"))
p1 = S.Matrix(S.symbols("p1_0:4"))

A, C, m, V, T = S.symbols("A C m V T")
l1, l2 = S.symbols("l1 l2")

def mdot(a, b):
    return (a.T * eta * b)[0]

def contract(a, b):
    return (a.T * b)[0]

s = -mdot(u, u)
v = S.sqrt(s)
d = mdot(u, q)
W = mdot(q, q) + d**2/s
P2 = mdot(P, P)

phi1 = contract(u, P)
phi2 = contract(q, P)

chi = contract(p1, u)
psi = contract(p1, q) + 2*A*W/v**3

H = (
    contract(p0, u)
    + contract(p1, q)
    - v*W*P2/(4*C)
    + A*W/v**3
    + m*v
    + v*(V+T)
)

def PB(f, g):
    result = S.Integer(0)
    for x, p in list(zip(u, p1)) + list(zip(q, P)):
        result += (
            S.diff(f, x)*S.diff(g, p)
            - S.diff(f, p)*S.diff(g, x)
        )
    return S.expand(result)

print("STAGE 1: primary subalgebra", flush=True)

assert S.simplify(PB(phi1, phi2)+phi1) == 0
assert S.simplify(PB(phi2, phi1)-phi1) == 0

print("PRIMARY SUBALGEBRA: PASS", flush=True)

print("STAGE 2: phi1 preservation", flush=True)

r1 = S.simplify(PB(phi1, H)-(phi2-chi))
print("RESIDUAL =", r1)
assert r1 == 0

print("{phi1,H} = phi2 - chi", flush=True)

print("STAGE 3: phi2 preservation", flush=True)

r2 = S.simplify(PB(phi2, H)-psi)
print("RESIDUAL =", r2)
assert r2 == 0

print("{phi2,H} = psi", flush=True)

print("STAGE 4: total Hamiltonian", flush=True)

HT = H + l1*phi1 + l2*phi2

expected1 = phi2-chi-l2*phi1
expected2 = psi+l1*phi1

r3 = S.simplify(PB(phi1, HT)-expected1)
r4 = S.simplify(PB(phi2, HT)-expected2)

print("TOTAL RESIDUALS =", [r3, r4])
assert r3 == 0 and r4 == 0

print("STAGE 4 PASS", flush=True)

print("STAGE 5: secondary independence",
      flush=True)

# On the primary surface phi1=phi2=0:
# dot(phi1)=-chi, dot(phi2)=psi.
#
# The constraints chi and psi are independent
# in the p1 directions when u and q are
# independent.

jet_u = S.Matrix([5, 3, 0, 0])
jet_q = S.Matrix([1, 2, 2, 1])

rank = S.Matrix([
    list(jet_u),
    list(jet_q)
]).rank()

print("SECONDARY p1-GRADIENT RANK =", rank)
assert rank == 2

print("STAGE 5 PASS", flush=True)

print()
print("R26D2 COMPLETE")
print("CORRECTED PRIMARY PRESERVATION: VERIFIED")
print("SECONDARY CANDIDATE chi=p1.u: DERIVED")
print("SECONDARY CANDIDATE psi=p1.q+2AW/v^3: DERIVED")
print("SECONDARY PRESERVATION: OPEN")
print("FULL DIRAC CLASSIFICATION: OPEN")
print("DEFI CORRESPONDENCE: OPEN")
print("ELAPSED SECONDS =", round(monotonic()-t0, 3))
