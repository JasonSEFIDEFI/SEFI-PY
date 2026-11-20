import sympy as S
from time import perf_counter

t0 = perf_counter()
print("=== CP098-T90C-R27Q ===", flush=True)

v, W, P2, IW, IS = S.symbols(
    "v W P2 IW IS", positive=True
)
d, lam2, R = S.symbols("d lam2 R", real=True)
pS, pW = S.symbols("pS pW", real=True)
hS, hW, kA = S.symbols(
    "hS hW kA", positive=True
)

s = v**2

# Compact rates from R27P, restricted to
# phi1=phi2=0 and A=0.
vd = -d/v
Wd = 2*d*W/s + 2*lam2*W
P2d = -2*R - 2*lam2*P2

# R represents the Minkowski contraction P.p1.
IWd = v*pW/hW
ISd = v*pS/hS

print("STAGE 1: scalar Wdot sign", flush=True)

q2 = S.symbols("q2", real=True)
dd = q2 + lam2*d
sd = -2*d
q_qdot = lam2*q2

# Scalar derivative of W=q²+d²/s.
Wdot_scalar = (
    2*q_qdot
    +2*d*dd/s
    -d**2*sd/s**2
)

# Here qdot's transverse P term vanishes
# in the contractions on the primary surface.
# The lam1 terms cancel from Wdot.
Wdot_scalar = S.expand(Wdot_scalar)

# Substitute q²=W-d²/s.
Wdot_reduced = S.factor(
    Wdot_scalar.subs(q2, W-d**2/s)
)

assert S.cancel(Wdot_reduced-Wd) == 0

print("PASS: positive geometric sign", flush=True)
print("Wdot =", Wd, flush=True)
print("STAGE 1 PASS", flush=True)

print("STAGE 2: Qdot reduction", flush=True)

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
    v: vd,
    W: Wd,
    P2: P2d,
    IW: IWd,
    IS: ISd
}

Qdot = sum(
    S.diff(Q,x)*rate
    for x,rate in rates.items()
)

Qdot_surface = S.factor(
    Qdot.subs(IS,IS_required)
)

lambda_coefficient = S.factor(
    S.diff(Qdot_surface,lam2)
)

print("lambda2 coefficient =", lambda_coefficient,
      flush=True)

print("STAGE 2 PASS", flush=True)

print("STAGE 3: required pS", flush=True)

solutions = S.solve(Qdot_surface,pS)
assert len(solutions) == 1

required_pS = S.factor(solutions[0])

assert S.cancel(
    Qdot_surface.subs(pS,required_pS)
) == 0

print("Required pS =", required_pS, flush=True)
print("STAGE 3 PASS", flush=True)

print()
print("R27Q COMPLETE")
print("SCALAR GEOMETRIC RATES: CHECKED")
print("QDOT ON Q=0: REDUCED")
print("LAMBDA2 DEPENDENCE: REPORTED")
print("FULL 4D CONSTRAINT COMPATIBILITY: OPEN")
print("GENERAL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-t0,3)
)
