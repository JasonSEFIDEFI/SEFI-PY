import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27N ===", flush=True)

IA, IS, IW = S.symbols("IA IS IW", real=True)
pA, pS = S.symbols("pA pS", real=True)
hA, hS = S.symbols("hA hS", positive=True)
kA, kS = S.symbols("kA kS", positive=True)
v, W = S.symbols("v W", positive=True)
P2 = S.symbols("P2", nonnegative=True)

A = (3 + 2*IA)*IS
C = (2 + IA + S.Rational(3,2)*IA**2)*IW

# Only the relevant identity-sector terms are needed.
H = (
    -v*W*P2/(4*C)
    +A*W/v**3
    +v*(kA*IA**2+kS*IS**2)/2
    +v*(pA**2/hA+pS**2/hS)/2
)

IAdot = S.diff(H, pA)
ISdot = S.diff(H, pS)
pAdot = -S.diff(H, IA)
pSdot = -S.diff(H, IS)

def verify(label, expr):
    residual = S.cancel(expr)
    if residual != 0:
        raise AssertionError(f"{label}: {residual}")
    print("PASS:", label, flush=True)

print("STAGE 1: identity-coordinate evolution", flush=True)

verify("IAdot", IAdot-v*pA/hA)
verify("ISdot", ISdot-v*pS/hS)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: fixed intersection", flush=True)

intersection = {
    IA: -S.Rational(3, 2),
    IS: 0,
    pA: 0,
    pS: 0
}

verify("IAdot at intersection",
       IAdot.subs(intersection))
verify("ISdot at intersection",
       ISdot.subs(intersection))

Q = S.factor(pAdot.subs(intersection))
verify(
    "Q identity",
    Q - (
        S.Rational(3,2)*v*kA
        +56*v*W*P2/(961*IW)
    )
)

verify(
    "pSdot at intersection",
    pSdot.subs(intersection)
)

print("Q =", Q, flush=True)
print("pSdot =", pSdot.subs(intersection),
      flush=True)

# At pA=0, the first derivative of v*pA/hA
# has no vdot contribution.
IA_ddot = S.factor(v*Q/hA)

print("IA second derivative =", IA_ddot,
      flush=True)

positive_example = {
    v: 2, W: 1, P2: 1,
    IW: 1, kA: 1, hA: 1
}
value = S.simplify(
    IA_ddot.subs(positive_example)
)
assert value > 0

print("Positive example =", value, flush=True)
print("STAGE 2 PASS", flush=True)

print("STAGE 3: IA=-3/2 branch", flush=True)

branch = {
    IA: -S.Rational(3,2),
    pA: 0
}

verify(
    "branch IAdot",
    IAdot.subs(branch)
)

branch_Q = S.factor(pAdot.subs(branch))

required_IS = (
    3*kA*v**4/(4*W)
    +28*P2*v**4/(961*IW)
)

verify(
    "branch Q vanishes on necessary surface",
    branch_Q.subs(IS, required_IS)
)

print("Necessary IS =", required_IS, flush=True)
print("Remaining requirement: preserve Q=0",
      flush=True)
print("STAGE 3 PASS", flush=True)

print()
print("R27N COMPLETE")
print("FIXED INTERSECTION: SECOND-ORDER OBSTRUCTION")
print("CONDITION: IW>0, kA>0, P2>=0")
print("IA BRANCH: Q=0 NECESSARY")
print("Q TANGENCY: OPEN")
print("B PRESERVATION: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3)
)
