import sympy as S
from time import monotonic

t0 = monotonic()
print("=== CP098-T90C-R26R ===", flush=True)

eta = S.diag(-1, 1, 1, 1)

u = S.Matrix(S.symbols("u0:4"))
q = S.Matrix(S.symbols("q0:4"))
P = S.Matrix(S.symbols("P0:4"))
p0 = S.Matrix(S.symbols("p0_0 p0_1 p0_2 p0_3"))
p1 = S.Matrix(S.symbols("p1_0 p1_1 p1_2 p1_3"))

def mdot(a, b):
    return (a.T * eta * b)[0]

def ordinary(a, b):
    return (a.T * b)[0]

def grad(expr, variables):
    return S.Matrix([
        S.diff(expr, x) for x in variables
    ])

def direction(expr, variables, vector):
    return ordinary(grad(expr, variables), vector)

def check(label, residual):
    result = S.cancel(residual)
    if result != 0:
        print("FAIL:", label, flush=True)
        print("RESIDUAL:", result, flush=True)
        raise AssertionError(label)
    print("PASS:", label, flush=True)

s = -mdot(u, u)
d = mdot(u, q)
W = mdot(q, q) + d**2/s

# Identity-sector variables.
I = S.Matrix(S.symbols("IA IS IW"))
pI = S.Matrix(S.symbols("pIA pIS pIW"))

mu1 = S.Function("mu1")
mu2 = S.Function("mu2")

A = mu1(I[0]) * I[1]
C = mu2(I[0]) * I[2]

h = S.symbols("hA hS hW", nonzero=True)
k = S.symbols("kA kS kW")
m = S.symbols("m")

V = sum(
    k[i] * I[i]**2 / 2
    for i in range(3)
)

T = sum(
    pI[i]**2 / (2*h[i])
    for i in range(3)
)

M = m + V + T

# v=sqrt(s), with s>0.
v = S.sqrt(s)

P2 = mdot(P, P)

D = -v*W*P2/(4*C)
F = A*W/v**3
G = v*M

Hc = (
    ordinary(p0, u)
    + ordinary(p1, q)
    + D + F + G
)

E = 2*F
psi = ordinary(q, p1) + E

phi1 = ordinary(u, P)
phi2 = ordinary(q, P)

print("STAGE 1: invariant derivatives", flush=True)

check("q.du(s)=-2d",
      direction(s, u, q)+2*d)

check("q.du(W)=2dW/s",
      direction(W, u, q)-2*d*W/s)

check("q.du(v)=-d/v",
      direction(v, u, q)+d/v)

check("q.du(vW)=dW/v",
      direction(v*W, u, q)-d*W/v)

check("q.du(E)=10AdW/v5",
      direction(E, u, q)-10*A*d*W/v**5)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: canonical decomposition", flush=True)

# Canonical Poisson bracket:
#
# {psi,Hc} =
#   (d_u psi).(d_p1 Hc)
# - (d_p1 psi).(d_u Hc)
# + (d_q psi).(d_P Hc)
# + (d_I psi).(d_pI Hc)
#
# psi has no dependence on P or p0.
# Hc has no dependence on y.
#
# The four contributions are computed
# independently below.

term_u_positive = ordinary(
    grad(psi, u), grad(Hc, p1)
)

term_u_negative = -ordinary(
    grad(psi, p1), grad(Hc, u)
)

term_qP = ordinary(
    grad(psi, q), grad(Hc, P)
)

term_identity = ordinary(
    grad(psi, I), grad(Hc, pI)
)

check(
    "positive u/p1 contribution",
    term_u_positive - 10*A*d*W/v**5
)

expected_negative = (
    -ordinary(p0, q)
    + d*W*P2/(4*C*v)
    - 5*A*d*W/v**5
    + d*M/v
)

check(
    "negative u/p1 contribution",
    term_u_negative - expected_negative
)

expected_qP = (
    -v*W*mdot(p1, P)/(2*C)
    -2*A*W/(C*v**2)
       *(phi2+d*phi1/s)
)

check(
    "q/P contribution",
    term_qP - expected_qP
)

expected_identity = (
    2*W/v**2
    *sum(
        S.diff(A, I[i])*pI[i]/h[i]
        for i in range(3)
    )
)

check(
    "identity-sector contribution",
    term_identity - expected_identity
)

print("STAGE 2 PASS", flush=True)

print("STAGE 3: assemble general B", flush=True)

B_direct = (
    term_u_positive
    + term_u_negative
    + term_qP
    + term_identity
)

B_expected = (
    -ordinary(p0, q)
    -v*W*mdot(p1, P)/(2*C)
    +d*W*P2/(4*C*v)
    +5*A*d*W/v**5
    +d*M/v
    -2*A*W/(C*v**2)
       *(phi2+d*phi1/s)
    +expected_identity
)

check(
    "GENERAL B DIRECT-FORMULA RESIDUAL",
    B_direct-B_expected
)

print("STAGE 3 PASS", flush=True)

print()
print("R26R COMPLETE")
print("FOUR CANONICAL CONTRIBUTIONS: PASS")
print("IDENTITY-SECTOR COUPLING: INCLUDED")
print("GENERAL B FORMULA: SYMBOLIC PASS")
print("NO NUMERICAL WITNESSES USED")
print("CONSTRAINT REGULARITY: OPEN")
print("A-ZERO BRANCH: OPEN")
print("GAUGE GENERATORS: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(monotonic()-t0, 3)
)
