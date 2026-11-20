import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27V ===", flush=True)

source = Path(__file__).with_name(
    "CP098_T90C_R27U.py"
)
assert source.is_file(), "R27U source missing"

# Reuse the previously verified symbolic definitions.
with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

z = S.symbols("z", real=True)
t = ns["t"]
point = dict(ns["point"])

p0 = ns["p0"]
p1 = ns["p1"]

point[p1[2]] = z
point[p0[1]] = -S.Rational(8,31)*z

def at(expr):
    return S.cancel(expr.subs(point))

print("STAGE 1: exact constraint preservation",
      flush=True)

constraints = [
    ("phi1",ns["phi1"]),
    ("phi2",ns["phi2"]),
    ("chi",ns["chi"]),
    ("psi",ns["psi"]),
    ("theta",ns["theta"]),
    ("B",ns["B"]),
    ("A",ns["A"]),
    ("J",ns["J"]),
    ("Q",ns["Q"]),
    ("pA",ns["pA"])
]

for name, expr in constraints:
    result = at(expr)
    assert result == 0, (name,result)
    print("PASS:",name)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: exact D(t,z)", flush=True)

pairs = (
    list(zip(ns["u"],ns["p1"]))
    + list(zip(ns["q"],ns["P"]))
    + [
        (ns["IA"],ns["pA"]),
        (ns["IS"],ns["pS"]),
        (ns["IW"],ns["pW"])
    ]
)

parts = []

for x,p in pairs:
    term = S.cancel(
        at(S.diff(ns["B"],x))
        *at(S.diff(ns["H"],p))
        -
        at(S.diff(ns["B"],p))
        *at(S.diff(ns["H"],x))
    )
    parts.append(term)

D = S.factor(sum(parts))

baseline = (
    3813248*t**2+593952197
)/S.Integer(14776336)

assert S.cancel(D.subs(z,0)-baseline) == 0

print("PASS: R27U recovered at z=0")
print("D(t,z) =",D)
print("STAGE 2 PASS", flush=True)

print("STAGE 3: search D=0 at t=0",
      flush=True)

D_z = S.factor(D.subs(t,0))
poly = S.Poly(D_z,z)

print("D(0,z) =",D_z)
print("degree in z =",poly.degree())

roots = S.solveset(
    S.Eq(D_z,0),z,domain=S.S.Reals
)

print("real roots =",roots)

for root in (
    list(roots) if isinstance(roots,S.FiniteSet)
    else []
):
    assert S.simplify(D_z.subs(z,root)) == 0

print("STAGE 3 PASS", flush=True)

print()
print("R27V COMPLETE")
print("TWO-PARAMETER FAMILY: VERIFIED")
print("D=0 CANDIDATES: REPORTED")
print("FULL BRANCH CLASSIFICATION: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(perf_counter()-start,3))
