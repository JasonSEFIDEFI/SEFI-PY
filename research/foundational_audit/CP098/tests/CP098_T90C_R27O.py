import sympy as S
from time import perf_counter

t0 = perf_counter()
print("=== CP098-T90C-R27O ===", flush=True)

v, W, P2, IW, IS = S.symbols(
    "v W P2 IW IS", real=True, nonzero=True
)
kA, hS, hW = S.symbols(
    "kA hS hW", positive=True
)
pS, pW = S.symbols("pS pW", real=True)
vd, Wd, P2d = S.symbols("vd Wd P2d", real=True)

Q = (
    S.Rational(3,2)*v*kA
    +56*v*W*P2/(961*IW)
    -2*IS*W/v**3
)

IS_required = (
    3*kA*v**4/(4*W)
    +28*P2*v**4/(961*IW)
)

print("STAGE 1: verify Q=0 surface", flush=True)

assert S.cancel(Q.subs(IS, IS_required)) == 0
print("PASS: Q=0 exactly", flush=True)

print("STAGE 2: derive total derivative", flush=True)

rates = {
    v: vd,
    W: Wd,
    P2: P2d,
    IW: v*pW/hW,
    IS: v*pS/hS
}

def total_derivative(F):
    return S.factor(sum(
        S.diff(F, x)*rate
        for x, rate in rates.items()
    ))

Qdot = total_derivative(Q)
restricted_Qdot = S.factor(
    Qdot.subs(IS, IS_required)
)

# On Q=0, differentiating the defining surface
# must agree with differentiating IS-IS_required.
surface_rate = S.factor(
    rates[IS] - total_derivative(IS_required)
)

assert S.cancel(
    restricted_Qdot
    +2*W*surface_rate/v**3
) == 0

print("PASS: equivalent tangency derivatives", flush=True)

print("STAGE 3: solve necessary pS condition",
      flush=True)

solutions = S.solve(surface_rate, pS)
assert len(solutions) == 1

pS_required = S.factor(solutions[0])

assert S.cancel(
    restricted_Qdot.subs(pS, pS_required)
) == 0

print("Required pS =", pS_required, flush=True)
print("PASS: Qdot=0 on restricted surface",
      flush=True)

print("STAGE 4: geometric-rate dependence",
      flush=True)

for rate in (vd, Wd, P2d):
    coefficient = S.factor(
        S.diff(pS_required, rate)
    )
    print(f"coefficient of {rate} =", coefficient,
          flush=True)

print("STAGE 4 PASS", flush=True)
print()
print("R27O COMPLETE")
print("Q=0 TANGENCY RELATION: DERIVED")
print("GEOMETRIC EVOLUTION RATES: UNRESOLVED")
print("FULL HAMILTONIAN MULTIPLIER COMPATIBILITY: OPEN")
print("B=0 PRESERVATION: OPEN")
print("GENERAL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(perf_counter()-t0,3))
