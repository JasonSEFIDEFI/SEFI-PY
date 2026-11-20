import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27M ===", flush=True)

IA, IS, IW = S.symbols("IA IS IW", real=True)
pA, pS, pW = S.symbols("pA pS pW", real=True)
hA, hS, hW = S.symbols("hA hS hW", positive=True)
kA, kS, kW = S.symbols("kA kS kW", positive=True)
v, W = S.symbols("v W", positive=True)
P2 = S.symbols("P2", nonnegative=True)

A = (3 + 2*IA)*IS
C = (2 + IA + S.Rational(3,2)*IA**2)*IW

H = (
    -v*W*P2/(4*C)
    +A*W/v**3
    +v*(
        kA*IA**2/2 + kS*IS**2/2 + kW*IW**2/2
        +pA**2/(2*hA)
        +pS**2/(2*hS)
        +pW**2/(2*hW)
    )
)

pAdot = -S.diff(H, IA)
pSdot = -S.diff(H, IS)

def check(label, expr):
    residual = S.cancel(expr)
    if residual != 0:
        raise AssertionError(
            f"{label}: residual = {residual}"
        )
    print("PASS:", label, flush=True)

print("STAGE 1: exact identity momentum equations",
      flush=True)

expected_pAdot = (
    -v*W*P2*S.diff(C,IA)/(4*C**2)
    -2*IS*W/v**3
    -v*kA*IA
)

expected_pSdot = (
    -(3+2*IA)*W/v**3
    -v*kS*IS
)

check("pAdot", pAdot-expected_pAdot)
check("pSdot", pSdot-expected_pSdot)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: IA=-3/2 component", flush=True)

branch = {IA:-S.Rational(3,2)}

branch_pAdot = S.factor(pAdot.subs(branch))

expected = (
    S.Rational(3,2)*v*kA
    +56*v*W*P2/(961*IW)
    -2*IS*W/v**3
)

check("branch pAdot", branch_pAdot-expected)

solutions = S.solve(branch_pAdot, IS)

if len(solutions) != 1:
    raise AssertionError("Unexpected IS solutions")

IS_required = S.factor(solutions[0])

check(
    "necessary IS condition",
    IS_required-(
        3*kA*v**4/(4*W)
        +28*P2*v**4/(961*IW)
    )
)

print("pAdot =", branch_pAdot, flush=True)
print("Required IS =", IS_required, flush=True)
print("STAGE 2 PASS", flush=True)

print("STAGE 3: intersection", flush=True)

intersection = {
    IA:-S.Rational(3,2),
    IS:0
}

Q = S.factor(pAdot.subs(intersection))

expected_Q = (
    S.Rational(3,2)*v*kA
    +56*v*W*P2/(961*IW)
)

check("intersection Q", Q-expected_Q)

print("Q =", Q, flush=True)

# This positivity result requires IW>0,
# in addition to the symbol assumptions.
positive_witness = {
    v:2, W:1, IW:1, P2:1, kA:1
}

q_value = S.simplify(Q.subs(positive_witness))
assert q_value > 0

print("Positive-parameter witness Q =", q_value,
      flush=True)
print(
    "For IW>0, kA>0, P2>=0, v,W>0: Q>0",
    flush=True
)

print("STAGE 3 PASS", flush=True)

print()
print("R27M COMPLETE")
print("IDENTITY MOMENTUM EQUATIONS: VERIFIED")
print("IA BRANCH NECESSARY CONDITION: DERIVED")
print("INTERSECTION Q POSITIVITY: CONDITIONAL")
print("HIGHER-ORDER TANGENCY: OPEN")
print("B=0 PRESERVATION: OPEN")
print("GENERAL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start,3)
)
