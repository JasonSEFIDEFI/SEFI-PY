import sympy as S
from time import perf_counter

t0 = perf_counter()
print("=== CP098-T90C-R27P ===", flush=True)

# Independent four-dimensional canonical construction.
eta = S.diag(-1, 1, 1, 1)

u = S.Matrix(S.symbols("u0:4", real=True))
q = S.Matrix(S.symbols("q0:4", real=True))
P = S.Matrix(S.symbols("P0:4", real=True))
p1 = S.Matrix(S.symbols("p10:14", real=True))

lam1, lam2 = S.symbols("lam1 lam2", real=True)
C, A, M = S.symbols("C A M", real=True, nonzero=True)

def dot(x, y):
    return (x.T * eta * y)[0]

s = -dot(u, u)
d = dot(u, q)
v = S.sqrt(s)
W = dot(q, q) + d**2/s
P2 = dot(P, P)

phi1 = (u.T*P)[0]
phi2 = (q.T*P)[0]

Hc = (
    (p1.T*q)[0]
    -v*W*P2/(4*C)
    +A*W/v**3
    +v*M
)

HT = Hc + lam1*phi1 + lam2*phi2

print("STAGE 1: canonical evolution", flush=True)

udot = S.Matrix([
    S.diff(HT, p1[i]) for i in range(4)
])

qdot = S.Matrix([
    S.diff(HT, P[i]) for i in range(4)
])

Pdot = S.Matrix([
    -S.diff(HT, q[i]) for i in range(4)
])

expected_udot = q
expected_qdot = (
    -v*W*eta*P/(2*C)
    +lam1*u+lam2*q
)

for i in range(4):
    assert S.simplify(
        udot[i]-expected_udot[i]
    ) == 0
    assert S.simplify(
        qdot[i]-expected_qdot[i]
    ) == 0

print("PASS: udot and qdot", flush=True)

print("STAGE 2: derive vdot", flush=True)

vdot = -d/v

# sdot=-2d follows from udot=q.
assert S.simplify(
    sum(S.diff(v,u[i])*udot[i]
        for i in range(4))-vdot
) == 0

print("vdot =", vdot, flush=True)
print("STAGE 2 PASS", flush=True)

print("STAGE 3: derive Wdot", flush=True)

Wdot_direct = sum(
    S.diff(W,u[i])*udot[i]
    +S.diff(W,q[i])*qdot[i]
    for i in range(4)
)

q2 = dot(q,q)
q_qdot = dot(q,qdot)
u_qdot = dot(u,qdot)

sdot = -2*d
ddot = q2+u_qdot

Wdot_scalar = (
    2*q_qdot
    +2*d*ddot/s
    -d**2*sdot/s**2
)

assert S.simplify(
    Wdot_direct-Wdot_scalar
) == 0

print("PASS: Wdot direct/scalar agreement",
      flush=True)

print("STAGE 4: derive P²dot", flush=True)

P2dot_direct = 2*dot(P,Pdot)

# Independently differentiate P² along the
# canonical P evolution.
P2dot_chain = sum(
    S.diff(P2,P[i])*Pdot[i]
    for i in range(4)
)

assert S.simplify(
    P2dot_direct-P2dot_chain
) == 0

print("PASS: P²dot agreement", flush=True)

print("STAGE 5: compact scalar formulas",
      flush=True)

# Since phi1=u.P and phi2=q.P:
q_qdot_expected = (
    -v*W*phi2/(2*C)
    +lam1*d+lam2*q2
)

u_qdot_expected = (
    -v*W*phi1/(2*C)
    -lam1*s+lam2*d
)

assert S.simplify(
    q_qdot-q_qdot_expected
) == 0
assert S.simplify(
    u_qdot-u_qdot_expected
) == 0

Wdot_compact = S.factor(
    Wdot_scalar.subs({
        q_qdot: q_qdot_expected,
        u_qdot: u_qdot_expected
    })
)

# Use direct substitution of the scalar
# contractions to avoid pattern ambiguity.
Wdot_compact = S.factor(
    2*q_qdot_expected
    +2*d*(q2+u_qdot_expected)/s
    +2*d**3/s**2
)

assert S.simplify(
    Wdot_direct-Wdot_compact
) == 0

print("PASS: compact Wdot identity",
      flush=True)

# Pdot can be expressed through dW/dq.
gradWq = S.Matrix([
    S.diff(W,q[i]) for i in range(4)
])

Pdot_expected = (
    -p1
    +(v*P2/(4*C)-A/v**3)*gradWq
    -lam2*P
)

for i in range(4):
    assert S.simplify(
        Pdot[i]-Pdot_expected[i]
    ) == 0

P2dot_compact = S.factor(
    2*dot(P,Pdot_expected)
)

assert S.simplify(
    P2dot_direct-P2dot_compact
) == 0

print("PASS: compact P²dot identity",
      flush=True)

print("STAGE 5 PASS", flush=True)

print()
print("R27P COMPLETE")
print("vdot =", vdot)
print("Wdot =", Wdot_compact)
print("P²dot =", P2dot_compact)
print("GEOMETRIC EVOLUTION RATES: DERIVED")
print("QDOT MULTIPLIER COMPATIBILITY: OPEN")
print("B PRESERVATION COMPATIBILITY: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-t0, 3)
)
