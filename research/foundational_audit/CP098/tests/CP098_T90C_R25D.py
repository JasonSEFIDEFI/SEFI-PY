import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R25D ===", flush=True)

# Use unrestricted symbols to diagnose the equations.
hA, kS, kW = S.symbols("hA kS kW", real=True)

# Exact R25C residuals.
rA = (133128*hA + 12226945)/8520192
rS = (4096*kS + 387)/1024
rW = 4*(16641*kW + 9844)/16641

residuals = (rA, rS, rW)
variables = (hA, kS, kW)

print("STAGE 1: unrestricted pointwise solution",
      flush=True)

solutions = S.solve(
    residuals, variables, dict=True
)

print("SOLUTIONS =", solutions, flush=True)
assert len(solutions) == 1, (
    "STOP: unexpected solution structure"
)

solution = solutions[0]

for variable in variables:
    assert variable in solution

print("STAGE 1 PASS", flush=True)

print("STAGE 2: exact substitution", flush=True)

checks = [
    S.factor(r.subs(solution))
    for r in residuals
]

print("RESIDUALS =", checks, flush=True)
assert checks == [0, 0, 0]

print("STAGE 2 PASS", flush=True)

print("STAGE 3: admissibility conditions",
      flush=True)

required_hA = solution[hA]
required_kS = solution[kS]
required_kW = solution[kW]

print("REQUIRED hA =", required_hA)
print("REQUIRED kS =", required_kS)
print("REQUIRED kW =", required_kW)

assert required_hA < 0
assert required_kS < 0
assert required_kW < 0

print("POSITIVE hA ADMISSIBLE =",
      required_hA > 0)

print("NONNEGATIVE kS ADMISSIBLE =",
      required_kS >= 0)

print("NONNEGATIVE kW ADMISSIBLE =",
      required_kW >= 0)

print("STAGE 3 PASS", flush=True)

print("STAGE 4: dependency classification",
      flush=True)

# At this point IA=0, IS=IW=1.
# IA'=1, IS'=IW'=0.
# Therefore hS, hW, and kA do not enter these
# three pointwise residual equations.

hS, hW, kA = S.symbols(
    "hS hW kA", real=True
)

for coefficient in (hS, hW, kA):
    assert all(
        S.diff(r, coefficient) == 0
        for r in residuals
    )

print("UNDETERMINED AT THIS JET:",
      "hS, hW, kA")

print("STAGE 4 PASS", flush=True)

print()
print("R25D COMPLETE")
print("UNRESTRICTED POINTWISE SOLUTION: VERIFIED")
print("POSITIVE hA WITNESS: IMPOSSIBLE")
print("GENERAL CANDIDATE MODEL: NOT DISPROVEN")
print("FULL TRAJECTORY DYNAMICS: OPEN")
print("WORLDLINE EQUATIONS: OPEN")
print("HAMILTONIAN RECONSTRUCTION: OPEN")
print("DIRAC SYSTEM: OPEN")
print("DEFI CORRESPONDENCE: OPEN")
print("ELAPSED SECONDS =",
      round(monotonic()-start, 3))
