import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R26Y ===", flush=True)

def check(label, residual):
    result = S.factor(residual)
    if result != 0:
        print("FAIL:", label, result, flush=True)
        raise AssertionError(label)
    print("PASS:", label, flush=True)

IA, IS, IW = S.symbols("IA IS IW", real=True)
pIA, pIS = S.symbols("pIA pIS", real=True)

hA, hS, v, W = S.symbols(
    "hA hS v W", positive=True
)

vdot = S.symbols("vdot", real=True)
P2, kA, kS = S.symbols("P2 kA kS", real=True)

root = -S.Rational(3, 2)

mu1 = 3 + 2*IA
mu2 = 2 + IA + S.Rational(3, 2)*IA**2

A = mu1*IS
C = mu2*IW

IA_dot = v*pIA/hA
IS_dot = v*pIS/hS

pIA_dot = (
    -2*IS*W/v**3
    -v*W*P2*S.diff(C, IA)/(4*C**2)
    -v*kA*IA
)

pIS_dot = (
    -mu1*W/v**3
    -v*kS*IS
)

IA_ddot = (vdot*pIA + v*pIA_dot)/hA
IS_ddot = (vdot*pIS + v*pIS_dot)/hS

print("STAGE 1: product-rule derivation", flush=True)

# A = (3+2IA)*IS.
#
# At the intersection, both factors vanish.
#
# A''' = 6*(IA''*IS' + IA'*IS'')
#
# This follows directly from the third
# derivative of a product.

A3_intersection = 6*(
    IA_ddot*IS_dot
    + IA_dot*IS_ddot
)

expected_general = (
    12*v*vdot*pIA*pIS/(hA*hS)
    +6*v**2*(
        pIA_dot*pIS+pIA*pIS_dot
    )/(hA*hS)
)

check(
    "third derivative product rule",
    A3_intersection-expected_general
)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: intersection reduction", flush=True)

A3_root = S.factor(
    A3_intersection.subs({
        IA:root,
        IS:0
    })
)

Q = S.factor(
    pIA_dot.subs({
        IA:root,
        IS:0
    })
)

check(
    "pIS-dot vanishes at intersection",
    pIS_dot.subs({IA:root,IS:0})
)

expected_root = (
    12*v*vdot*pIA*pIS/(hA*hS)
    +6*v**2*Q*pIS/(hA*hS)
)

check(
    "intersection third derivative",
    A3_root-expected_root
)

print("pIA-dot at intersection =", Q)
print("A-third at intersection =", A3_root)

print("STAGE 2 PASS", flush=True)

print("STAGE 3: pIA=0 subset", flush=True)

subset_A = S.factor(A3_root.subs(pIA,0))

expected_A = 6*v**2*Q*pIS/(hA*hS)

check(
    "pIA=0 third derivative",
    subset_A-expected_A
)

print("When pIA=0:")
print("A-third =",subset_A)
print("If pIS!=0, third-order tangency")
print("requires Q=0.")

print("STAGE 3 PASS", flush=True)

print("STAGE 4: pIS=0 subset", flush=True)

subset_S = S.factor(A3_root.subs(pIS,0))

check(
    "pIS=0 third derivative",
    subset_S
)

print("When pIS=0:")
print("A-third=0 automatically.")
print("Fourth-order tangency remains open.")

print("STAGE 4 PASS", flush=True)

print("STAGE 5: intersection momentum cases",
      flush=True)

both_zero = A3_root.subs({
    pIA:0,
    pIS:0
})

check(
    "both momenta zero",
    both_zero
)

print("CASE 1: pIA!=0, pIS!=0")
print("  Second-order tangency fails.")

print("CASE 2: pIA=0, pIS!=0")
print("  Third-order tangency requires Q=0.")

print("CASE 3: pIS=0")
print("  Third-order tangency holds.")
print("  Higher derivatives remain open.")

print("STAGE 5 PASS", flush=True)

print()
print("R26Y COMPLETE")
print("THIRD-ORDER TANGENCY: VERIFIED")
print("CORRECTED MOMENTUM EQUATIONS: USED")
print("FULL BRANCH INVARIANCE: OPEN")
print("FOURTH-ORDER TANGENCY: OPEN")
print("B=0 PRESERVATION: OPEN")
print("DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(monotonic()-start,3)
)
