import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R25F ===", flush=True)

# Exact quantities established in R25C/R25E.
v = S.Integer(4)
vp = -S.Rational(1, 4)

sourceA = -S.Rational(12226945, 8520192)
sourceW = -S.Rational(39376, 16641)

# Positive kinetic coefficients, with kW initially
# unrestricted so its required value can be derived.
hA, hW = S.symbols("hA hW", positive=True)
kW = S.symbols("kW", real=True)

print("STAGE 1: general identity accelerations",
      flush=True)

# At this jet:
# IA=0, IA'=1, IW=1, IW'=0.
#
# kA drops out because IA=0.
accA = S.factor(vp/v + v*sourceA/hA)
accW = S.factor(v*(sourceW-v*kW)/hW)

print("IA'' =", accA)
print("IW'' =", accW)

# Independently recover the R25E unit-coefficient
# accelerations.
unitA = S.factor(accA.subs(hA, 1))
unitW = S.factor(accW.subs({hW: 1, kW: 1}))

assert unitA == -S.Rational(12360073, 2130048)
assert unitW == -S.Rational(423760, 16641)

print("R25E ACCELERATION CROSS-CHECK: PASS",
      flush=True)

print("STAGE 2: C-second-derivative chain rule",
      flush=True)

lam = S.symbols("lam", real=True)
alpha, omega = S.symbols("alpha omega", real=True)

IA = lam + alpha*lam**2/2
IW = 1 + omega*lam**2/2

C = (
    2 + IA + S.Rational(3, 2)*IA**2
) * IW

C2 = S.factor(S.diff(C, lam, 2).subs(lam, 0))

assert S.expand(C2 - (3 + alpha + 2*omega)) == 0

print("GENERAL C'' =", C2)
print("TARGET CONDITION: IA'' + 2*IW'' = 0")
print("STAGE 2 PASS", flush=True)

print("STAGE 3: nonnegative-potential obstruction",
      flush=True)

# For hA>0:
# vp/v is negative, and sourceA is negative.
# Therefore accA<0.
#
# For hW>0 and kW>=0:
# sourceW-v*kW is negative.
# Therefore accW<0.
#
# Both inequalities are exact.

assert vp/v < 0
assert sourceA < 0
assert sourceW < 0
assert v > 0

print("IF hA>0: IA'' < 0")
print("IF hW>0 AND kW>=0: IW'' < 0")
print("THEREFORE IA'' + 2*IW'' < 0")
print("TARGET C''=3: IMPOSSIBLE UNDER THESE CONDITIONS")

print("STAGE 3 PASS", flush=True)

print("STAGE 4: unrestricted kW requirement",
      flush=True)

condition = S.factor(accA + 2*accW)
solutions = S.solve(condition, kW)

assert len(solutions) == 1

required_kW = S.factor(solutions[0])

print("REQUIRED kW =", required_kW)

# Verify algebraically for arbitrary positive hA,hW.
check = S.factor(
    condition.subs(kW, required_kW)
)

assert check == 0
print("SYMBOLIC RESIDUAL =", check)

# Exact example: positive unit kinetic coefficients.
example_kW = S.factor(
    required_kW.subs({hA: 1, hW: 1})
)

print("WITH hA=hW=1:")
print("REQUIRED kW =", example_kW)

assert example_kW < 0

example_A = unitA
example_W = S.factor(
    accW.subs({hW: 1, kW: example_kW})
)

example_C2 = S.factor(
    3 + example_A + 2*example_W
)

print("IA'' =", example_A)
print("IW'' =", example_W)
print("RESULTING C'' =", example_C2)

assert example_C2 == 3

print("STAGE 4 PASS", flush=True)

print()
print("R25F COMPLETE")
print("POSITIVE-KINETIC/NONNEGATIVE-kW OBSTRUCTION: VERIFIED")
print("UNRESTRICTED-kW POINTWISE WITNESS: VERIFIED")
print("FULL TRAJECTORY: NOT ESTABLISHED")
print("WORLDLINE EQUATIONS: NOT SOLVED")
print("FULL HAMILTONIAN: OPEN")
print("DIRAC SYSTEM: OPEN")
print("DEFI CORRESPONDENCE: OPEN")
print("ELAPSED SECONDS =", round(monotonic()-start, 3))
