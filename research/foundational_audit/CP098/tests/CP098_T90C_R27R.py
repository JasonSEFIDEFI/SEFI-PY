import sympy as S
from time import perf_counter

t0 = perf_counter()
print("=== CP098-T90C-R27R ===", flush=True)

v, W, IW = S.symbols("v W IW", positive=True)
P2 = S.symbols("P2", nonnegative=True)
kA, hS, hW = S.symbols("kA hS hW", positive=True)
d, R, pS, pW = S.symbols(
    "d R pS pW", real=True
)
lam2, D = S.symbols("lam2 D", real=True)

IS = S.symbols("IS", real=True)

Q = (
    S.Rational(3,2)*v*kA
    +56*v*W*P2/(961*IW)
    -2*IS*W/v**3
)

IS_required = (
    3*kA*v**4/(4*W)
    +28*P2*v**4/(961*IW)
)

rates = {
    v: -d/v,
    W: 2*W*d/v**2 + 2*W*lam2,
    P2: -2*R - 2*lam2*P2,
    IW: v*pW/hW,
    IS: v*pS/hS
}

def exact_zero(expr):
    return S.cancel(expr) == 0

print("STAGE 1: construct restricted Qdot",
      flush=True)

assert exact_zero(Q.subs(IS, IS_required))

Qdot = sum(
    S.diff(Q, x)*rate
    for x, rate in rates.items()
)

Qdot_branch = S.factor(
    Qdot.subs(IS, IS_required)
)

L = S.factor(S.diff(Qdot_branch, lam2))
F = S.factor(Qdot_branch.subs(lam2, 0))

assert exact_zero(Qdot_branch-(F+lam2*L))

expected_L = (
    -v*(2883*IW*kA+112*P2*W)/(961*IW)
)

assert exact_zero(L-expected_L)

print("PASS: Qdot = F + lambda2*L")
print("L =", L)
print("STAGE 1 PASS")

print("STAGE 2: multiplier solvability",
      flush=True)

lam2_required = S.factor(-F/L)

assert exact_zero(
    Qdot_branch.subs(lam2, lam2_required)
)

print("PASS: lambda2 = -F/L")
print("F =", F)
print("lambda2 =", lam2_required)
print(
    "Under IW,kA,v,W>0 and P2>=0: L<0"
)
print("STAGE 2 PASS")

print("STAGE 3: independent B preservation",
      flush=True)

# Previously derived weak identity:
# Bdot = D + lambda2*(2WJ/v^2).
# On IA=-3/2 and pA=0, J=0.

J = S.symbols("J", real=True)

Bdot_weak = D + lam2*2*W*J/v**2

Bdot_branch = S.simplify(
    Bdot_weak.subs(J, 0)
)

assert exact_zero(Bdot_branch-D)
assert exact_zero(
    S.diff(Bdot_branch, lam2)
)

print("PASS: Bdot = D when J=0")
print("PASS: lambda2 cannot cancel D")
print("STAGE 3 PASS")

print()
print("R27R COMPLETE")
print("QDOT MULTIPLIER SOLVABILITY: VERIFIED")
print("B PRESERVATION REQUIRES D=0")
print("D ON Q=0 BRANCH: NOT YET EVALUATED")
print("MULTIPLIER COMPATIBILITY: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-t0, 3)
)
