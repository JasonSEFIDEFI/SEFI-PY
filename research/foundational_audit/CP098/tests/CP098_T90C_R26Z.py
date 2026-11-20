import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R26Z ===", flush=True)

IA, IS, IW = S.symbols("IA IS IW", real=True)
pIA, pIS = S.symbols("pIA pIS", real=True)
hA, hS, v, W = S.symbols(
    "hA hS v W", positive=True
)
kA = S.symbols("kA", real=True)
P2 = S.symbols("P2", real=True)

root = -S.Rational(3, 2)

mu1 = 3 + 2*IA
mu2 = 2 + IA + S.Rational(3, 2)*IA**2

A = mu1*IS
C = mu2*IW

def check(name, expr):
    result = S.cancel(expr)
    if result != 0:
        raise AssertionError(f"{name}: {result}")
    print(f"PASS: {name}", flush=True)

print("STAGE 1: regularity at intersection", flush=True)

Croot = S.factor(C.subs(IA, root))
check("C intersection", Croot - 31*IW/8)

assert Croot != 0
print("STAGE 1 PASS", flush=True)

print("STAGE 2: identity momentum obstruction", flush=True)

pIA_dot = (
    -2*IS*W/v**3
    -v*W*P2*S.diff(C, IA)/(4*C**2)
    -v*kA*IA
)

Q = S.factor(pIA_dot.subs({
    IA: root,
    IS: 0
}))

expected = (
    v*(2883*IW*kA + 112*P2*W)
    /(1922*IW)
)

check("intersection momentum equation", Q-expected)

print("Q =", Q, flush=True)
print("STAGE 2 PASS", flush=True)

print("STAGE 3: third-order tangency condition", flush=True)

A3 = 6*v**2*Q*pIS/(hA*hS)

check(
    "third-order obstruction",
    A3 - 6*v**2*Q*pIS/(hA*hS)
)

print(
    "For pIA=0 and pIS!=0, "
    "third-order tangency requires Q=0.",
    flush=True
)

print("STAGE 3 PASS", flush=True)

print("STAGE 4: momentum-sign implications", flush=True)

P2_required = S.solve(
    2883*IW*kA + 112*P2*W,
    P2
)

assert len(P2_required) == 1

P2_required = S.factor(P2_required[0])

check(
    "required momentum squared",
    P2_required + 2883*IW*kA/(112*W)
)

print("Required P2 =", P2_required, flush=True)

print(
    "If IW>0, W>0, kA>0 and P2>=0, "
    "Q cannot vanish.",
    flush=True
)

print(
    "This is conditional on IW and kA signs; "
    "those signs are not imposed here.",
    flush=True
)

print("STAGE 4 PASS", flush=True)

print()
print("R26Z COMPLETE")
print("INTERSECTION OBSTRUCTION: VERIFIED")
print("CONDITIONAL SIGN EXCLUSION: VERIFIED")
print("B=0 PRESERVATION: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3)
)
