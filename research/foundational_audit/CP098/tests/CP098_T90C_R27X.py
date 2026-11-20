import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27X ===", flush=True)

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

# Compute derivative expressions once.
print("STAGE 1: canonical derivative preparation",
      flush=True)

derivatives = []
for x, p in pairs:
    derivatives.append((
        S.diff(ns["B"], x),
        S.diff(ns["B"], p),
        S.diff(ns["H"], x),
        S.diff(ns["H"], p)
    ))

print("PASS: 11 canonical pairs prepared")
print("STAGE 1 PASS", flush=True)

def witness(rho):
    rho = S.Rational(rho)
    point = dict(ns["point"])

    # Fix the inherited R27U free parameter t=0.
    point[p1[3]] = S.Integer(0)

    point[P[2]] = rho

    P2 = rho**2
    IS = 12 + S.Rational(448, 961)*P2

    point[ns["IS"]] = IS

    IA = -S.Rational(3, 2)
    IW = S.Integer(1)
    C = S.Rational(31, 8)
    v = S.Integer(2)
    W = S.Integer(1)

    M = 1 + (IA**2 + IS**2 + IW**2)/2

    # theta=0 determines p0[0].
    point[p0[0]] = (
        v*W*P2/(8*C) - v*M/2
    )

    return point

def evaluate(expr, point):
    return S.cancel(expr.subs(point))

print("STAGE 2: exact witness audit and D",
      flush=True)

results = {}

for rho in (0, 1, 2):
    print("rho =", rho, flush=True)
    point = witness(rho)

    for name, expr in constraints:
        residual = evaluate(expr, point)
        assert residual == 0, (
            rho, name, residual
        )

    print("PASS: all ten constraints", flush=True)

    terms = []
    for Bx, Bp, Hx, Hp in derivatives:
        term = S.cancel(
            evaluate(Bx, point)*evaluate(Hp, point)
            -evaluate(Bp, point)*evaluate(Hx, point)
        )
        terms.append(term)

    D = S.cancel(sum(terms))
    results[rho] = D

    print("P² =", rho**2)
    print("IS =", point[ns["IS"]])
    print("D =", D)
    print("sign =", S.sign(D), flush=True)

baseline = S.Rational(
    593952197, 14776336
)
assert results[1] == baseline

print("PASS: R27T baseline recovered")
print("STAGE 2 PASS", flush=True)

print("STAGE 3: interpretation", flush=True)

for rho, D in results.items():
    if D == 0:
        print(
            f"rho={rho}: survives D=0 condition"
        )
    else:
        print(
            f"rho={rho}: excluded by D != 0"
        )

print("STAGE 3 PASS")

print()
print("R27X COMPLETE")
print("GEOMETRIC MOMENTUM WITNESSES: CHECKED")
print("GENERAL P² DEPENDENCE: OPEN")
print("GENERAL IA BRANCH: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3)
)
