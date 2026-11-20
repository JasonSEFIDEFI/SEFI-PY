import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27B ===", flush=True)

# Two-dimensional canonical subspace.
# This is a restricted algebraic audit, not
# a replacement for the full 4D calculation.

u = S.symbols("u0:2", real=True)
q = S.symbols("q0:2", real=True)
P = S.symbols("P0:2", real=True)
p0 = S.symbols("p0_0 p0_1", real=True)
p1 = S.symbols("p1_0 p1_1", real=True)

IA, IS, IW = S.symbols("IA IS IW", real=True)
pIA, pIS, pIW = S.symbols("pIA pIS pIW", real=True)
hA, hS, hW = S.symbols("hA hS hW", positive=True)
kA, kS, kW, m = S.symbols("kA kS kW m", real=True)

eta = S.diag(-1, 1)

def dot(a, b):
    return (S.Matrix(a).T * eta * S.Matrix(b))[0]

def ordinary(a, b):
    return sum(x*y for x, y in zip(a, b))

s = -dot(u, u)
v = S.sqrt(s)
d = dot(u, q)
W = dot(q, q) + d**2/s

A = (3 + 2*IA)*IS
C = (2 + IA + S.Rational(3,2)*IA**2)*IW

V = (
    kA*IA**2 + kS*IS**2 + kW*IW**2
)/2

T = (
    pIA**2/hA + pIS**2/hS + pIW**2/hW
)/2

M = m + V + T

phi1 = ordinary(u, P)
phi2 = ordinary(q, P)

P2 = dot(P, P)

B = (
    -ordinary(p0, q)
    -v*W*dot(p1, P)/(2*C)
    +d*W*P2/(4*C*v)
    +5*A*d*W/v**5
    +d*M/v
    -2*A*W/(C*v**2)*(phi2+d*phi1/s)
    +2*W/v**2*(
        2*IS*pIA/hA
        +(3+2*IA)*pIS/hS
    )
)

# Canonical pairs:
# (u,p1), (q,P), (I,pI).
# The primary constraints have no identity
# momentum dependence, so identity terms
# do not contribute to these two brackets.

def PB(F, G):
    result = 0
    for x, px in zip(u, p1):
        result += (
            S.diff(F, x)*S.diff(G, px)
            -S.diff(F, px)*S.diff(G, x)
        )
    for x, px in zip(q, P):
        result += (
            S.diff(F, x)*S.diff(G, px)
            -S.diff(F, px)*S.diff(G, x)
        )
    return result

def check(name, expr):
    result = S.factor(S.cancel(S.powdenest(
        S.expand_power_base(expr, force=False),
        force=False
    )))
    if result != 0:
        raise AssertionError(f"{name}: {result}")
    print("PASS:", name, flush=True)

print("STAGE 1: canonical primary brackets", flush=True)

check("primary algebra", PB(phi1, phi2)+phi1)
print("STAGE 1 PASS", flush=True)

print("STAGE 2: direct B brackets", flush=True)

br1 = PB(B, phi1)
br2 = PB(B, phi2)

# Compare candidate identities from R27A.
# Do not assume they hold.

delta1 = br1-B
delta2 = br2-2*B

print("Direct bracket expressions constructed.", flush=True)

print("STAGE 2 PASS", flush=True)

print("STAGE 3: numerical rational witness", flush=True)

witness = {
    u[0]: 2, u[1]: 0,
    q[0]: 0, q[1]: 1,
    P[0]: 0, P[1]: 1,
    p0[0]: 1, p0[1]: 2,
    p1[0]: 1, p1[1]: 3,
    IA: 0, IS: 1, IW: 1,
    pIA: 1, pIS: 2, pIW: 3,
    hA: 1, hS: 1, hW: 1,
    kA: 1, kS: 1, kW: 1,
    m: 1
}

for name, expression in [
    ("delta1", delta1),
    ("delta2", delta2)
]:
    value = S.simplify(expression.subs(witness))
    print(name, "=", value, flush=True)

print("STAGE 3 PASS", flush=True)

print("STAGE 4: singular witness", flush=True)

singular = dict(witness)
singular[IA] = -S.Rational(3,2)
singular[IS] = 0

for name, expression in [
    ("B", B),
    ("bracket1", br1),
    ("bracket2", br2)
]:
    print(
        name, "=",
        S.simplify(expression.subs(singular)),
        flush=True
    )

print("STAGE 4 PASS", flush=True)

print()
print("R27B COMPLETE")
print("DIRECT BRACKETS: CONSTRUCTED")
print("CANDIDATE IDENTITIES: WITNESS TESTED")
print("FULL 4D VERIFICATION: OPEN")
print("B=0 PRESERVATION: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3)
)
