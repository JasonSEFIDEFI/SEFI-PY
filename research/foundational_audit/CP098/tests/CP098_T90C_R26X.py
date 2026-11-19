import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R26X ===", flush=True)

IA, IS = S.symbols("IA IS", real=True)
pIA, pIS = S.symbols("pIA pIS", real=True)

hA, hS, v, W = S.symbols(
    "hA hS v W", positive=True
)

vdot = S.symbols("vdot", real=True)
P2, C, CIA, kA, kS = S.symbols(
    "P2 C CIA kA kS", real=True
)

def check(label, residual):
    r = S.factor(residual)
    if r != 0:
        print("FAIL:", label, r, flush=True)
        raise AssertionError(label)
    print("PASS:", label, flush=True)

A = (3 + 2*IA)*IS

IA_dot = v*pIA/hA
IS_dot = v*pIS/hS

pIA_dot = (
    -2*IS*W/v**3
    -v*W*P2*CIA/(4*C**2)
    -v*kA*IA
)

pIS_dot = (
    -(3+2*IA)*W/v**3
    -v*kS*IS
)

print("STAGE 1: first derivative", flush=True)

A_dot = (
    S.diff(A, IA)*IA_dot
    + S.diff(A, IS)*IS_dot
)

expected = v*(
    2*IS*pIA/hA
    + (3+2*IA)*pIS/hS
)

check("A-dot", A_dot-expected)
print("STAGE 1 PASS", flush=True)

print("STAGE 2: second derivative", flush=True)

# Full total derivative, including vdot.
A_ddot = (
    S.diff(A_dot, IA)*IA_dot
    + S.diff(A_dot, IS)*IS_dot
    + S.diff(A_dot, pIA)*pIA_dot
    + S.diff(A_dot, pIS)*pIS_dot
    + S.diff(A_dot, v)*vdot
)

root = -S.Rational(3, 2)

intersection = S.factor(
    A_ddot.subs({IA:root, IS:0})
)

expected_intersection = (
    4*v**2*pIA*pIS/(hA*hS)
)

check(
    "intersection second derivative",
    intersection-expected_intersection
)

print("A-ddot at intersection =", intersection)
print("STAGE 2 PASS", flush=True)

print("STAGE 3: exact obstruction", flush=True)

check(
    "pIA=0 makes A-ddot vanish",
    intersection.subs(pIA, 0)
)

check(
    "pIS=0 makes A-ddot vanish",
    intersection.subs(pIS, 0)
)

print("Necessary second-order condition:")
print("pIA*pIS = 0")
print("STAGE 3 PASS", flush=True)

print("STAGE 4: interpretation", flush=True)

print("At IA=-3/2 and IS=0:")
print("  A-dot = 0 automatically")
print("  A-ddot != 0 if pIA*pIS != 0")
print("  Higher-order tangency remains open")
print("  Preservation of B=0 remains open")

print("STAGE 4 PASS", flush=True)

print()
print("R26X COMPLETE")
print("SECOND-ORDER TANGENCY: VERIFIED")
print("INTERSECTION MOMENTUM CONDITION: VERIFIED")
print("FULL SINGULAR-BRANCH INVARIANCE: OPEN")
print("B=0 PRESERVATION: OPEN")
print("GAUGE GENERATORS: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(monotonic()-start, 3)
)
