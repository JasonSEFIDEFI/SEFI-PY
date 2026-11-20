import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27Z ===", flush=True)

source = Path(__file__).with_name(
    "CP098_T90C_R27U.py"
)
assert source.is_file(), "R27U source missing"

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

u, q, P = ns["u"], ns["q"], ns["P"]
p0, p1 = ns["p0"], ns["p1"]

pairs = (
    list(zip(u, p1))
    + list(zip(q, P))
    + [
        (ns["IA"], ns["pA"]),
        (ns["IS"], ns["pS"]),
        (ns["IW"], ns["pW"])
    ]
)

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

print("STAGE 1: canonical derivative preparation",
      flush=True)

derivatives = [
    (
        S.diff(ns["B"], x),
        S.diff(ns["B"], p),
        S.diff(ns["H"], x),
        S.diff(ns["H"], p)
    )
    for x, p in pairs
]

print("PASS: 11 canonical pairs prepared")
print("STAGE 1 PASS", flush=True)

def make_point(w):
    w = S.Rational(w)
    assert w != 0

    point = dict(ns["point"])
    point[q[1]] = w
    point[P[2]] = S.Integer(1)
    point[p1[3]] = S.Integer(0)

    W0 = w**2
    P20 = S.Integer(1)
    IS0 = 12/W0 + S.Rational(448, 961)

    point[ns["IS"]] = IS0

    IA0 = -S.Rational(3, 2)
    C0 = S.Rational(31, 8)
    v0 = S.Integer(2)

    M0 = 1 + (IA0**2 + IS0**2 + 1)/2

    # theta=0, with p1.q=0 and A=0.
    point[p0[0]] = (
        v0*W0*P20/(8*C0) - v0*M0/2
    )

    return point

def at(expr, point):
    return S.cancel(expr.subs(point))

print("STAGE 2: exact witness evaluations",
      flush=True)

results = {}

for w in (
    S.Rational(1, 2),
    S.Integer(1),
    S.Integer(2)
):
    point = make_point(w)

    for name, expr in constraints:
        residual = at(expr, point)
        assert residual == 0, (w, name, residual)

    print("PASS: all ten constraints at w =", w,
          flush=True)

    parts = []

    for Bx, Bp, Hx, Hp in derivatives:
        term = S.cancel(
            at(Bx, point)*at(Hp, point)
            - at(Bp, point)*at(Hx, point)
        )
        parts.append(term)

    D = S.cancel(sum(parts))
    results[w] = D

    print("W =", w**2)
    print("IS =", point[ns["IS"]])
    print("geometric contribution =",
          S.cancel(sum(parts[:8])))
    print("identity contribution =",
          S.cancel(sum(parts[8:])))
    print("D =", D)
    print("sign =", S.sign(D), flush=True)

baseline = S.Rational(
    593952197, 14776336
)

assert results[S.Integer(1)] == baseline

print("PASS: R27T baseline reproduced")
print("STAGE 2 PASS", flush=True)

print("STAGE 3: consistency classification",
      flush=True)

for w, D in results.items():
    if D == 0:
        print("w =", w, ": survives D=0")
    else:
        print("w =", w, ": fails D=0")

print("STAGE 3 PASS", flush=True)

print()
print("R27Z COMPLETE")
print("VARIABLE-W WITNESSES: CHECKED")
print("GENERAL W DEPENDENCE: OPEN")
print("GENERAL IA BRANCH: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3)
)
