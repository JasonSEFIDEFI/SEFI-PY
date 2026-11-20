import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27W ===", flush=True)

source = Path(__file__).with_name(
    "CP098_T90C_R27U.py"
)
assert source.is_file(), "R27U source missing"

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

a, b = S.symbols("a b", real=True)
point = dict(ns["point"])

point[ns["pS"]] = a
point[ns["pW"]] = b

# M increases by (a^2+b^2)/2.
# At u=(2,0,0,0), theta=0 requires
# p0[0] to decrease by (a^2+b^2)/2.
point[ns["p0"][0]] -= (a**2+b**2)/2

def at(expr):
    return S.cancel(expr.subs(point))

print("STAGE 1: constraint audit", flush=True)

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

print("STAGE 2: canonical D(a,b)", flush=True)

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

baseline = S.Rational(
    593952197, 14776336
)

assert S.cancel(D.subs({a:0,b:0})-baseline) == 0

print("PASS: R27T baseline reproduced")
print("D(a,b) =", D)
print("STAGE 2 PASS", flush=True)

print("STAGE 3: polynomial classification",
      flush=True)

poly = S.Poly(D, a, b)

print("total degree =", poly.total_degree())
print("D(a,0) =", S.factor(D.subs(b,0)))
print("D(0,b) =", S.factor(D.subs(a,0)))

# A zero at a=b=0 is already excluded.
# Analyze possible nonzero roots separately.
roots_a = S.solveset(
    S.Eq(D.subs(b,0),0),
    a, domain=S.S.Reals
)
roots_b = S.solveset(
    S.Eq(D.subs(a,0),0),
    b, domain=S.S.Reals
)

print("real roots along b=0:", roots_a)
print("real roots along a=0:", roots_b)

print("STAGE 3 PASS", flush=True)

print()
print("R27W COMPLETE")
print("IDENTITY MOMENTUM FAMILY: CHECKED")
print("D=0 ROOTS ON AXES: REPORTED")
print("FULL TWO-DIMENSIONAL ROOT SET: OPEN")
print("GENERAL BRANCH CONSISTENCY: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(perf_counter()-start,3))
