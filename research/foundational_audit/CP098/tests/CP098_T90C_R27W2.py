import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27W2 ===", flush=True)

source = Path(__file__).with_name(
    "CP098_T90C_R27U.py"
)
assert source.is_file(), "R27U source missing"

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

a, b = S.symbols("a b", real=True)
t = ns["t"]

point = dict(ns["point"])
point[ns["pS"]] = a
point[ns["pW"]] = b

# Preserve theta=0 as identity kinetic energy changes.
point[ns["p0"][0]] -= (a**2+b**2)/2

def at(expr):
    return S.cancel(expr.subs(point))

print("STAGE 1: exact constraint audit", flush=True)

constraints = [
    ("phi1", ns["phi1"]),
    ("phi2", ns["phi2"]),
    ("chi", ns["chi"]),
    ("psi", ns["psi"]),
    ("theta", ns["theta"]),
    ("B", ns["B"]),
    ("A", ns["A"]),
    ("J", ns["J"]),
    ("Q", ns["Q"]),
    ("pA", ns["pA"])
]

for name, expr in constraints:
    residual = at(expr)
    assert residual == 0, (name, residual)
    print("PASS:", name)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: canonical bracket", flush=True)

pairs = (
    list(zip(ns["u"], ns["p1"]))
    + list(zip(ns["q"], ns["P"]))
    + [
        (ns["IA"], ns["pA"]),
        (ns["IS"], ns["pS"]),
        (ns["IW"], ns["pW"])
    ]
)

parts = []

for i, (x, p) in enumerate(pairs, 1):
    term = S.cancel(
        at(S.diff(ns["B"], x))
        * at(S.diff(ns["H"], p))
        -
        at(S.diff(ns["B"], p))
        * at(S.diff(ns["H"], x))
    )
    parts.append(term)
    print(f"pair {i:02d} = {term}", flush=True)

D = S.factor(sum(parts))

baseline = S.Rational(593952197, 14776336)

expected = (
    baseline
    + (a**2+b**2)/4
    + 8*t**2/31
)

assert S.cancel(D-expected) == 0
assert S.cancel(D.subs({a:0,b:0,t:0})-baseline) == 0

print("PASS: full three-parameter identity")
print("PASS: R27T baseline at a=b=t=0")
print("D(a,b,t) =", D)
print("STAGE 2 PASS", flush=True)

print("STAGE 3: positivity certificate", flush=True)

assert baseline > 0
assert S.cancel(D-baseline) == (
    a**2/4+b**2/4+8*t**2/31
)

print("D - baseline =", S.factor(D-baseline))
print("baseline =", baseline)
print("PASS: D >= baseline > 0 for all real a,b,t")
print("PASS: no real D=0 solutions in this family")
print("STAGE 3 PASS", flush=True)

print()
print("R27W2 COMPLETE")
print("THREE-PARAMETER CONSTRAINT FAMILY: VERIFIED")
print("D POSITIVITY: EXACT")
print("D=0 SOLUTIONS IN THIS FAMILY: NONE")
print("GENERAL IA BRANCH: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start,3)
)
