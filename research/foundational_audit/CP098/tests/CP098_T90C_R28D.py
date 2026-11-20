import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28D ===", flush=True)

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

print("STAGE 1: canonical derivatives", flush=True)

derivatives = [
    (
        S.diff(ns["B"], x),
        S.diff(ns["B"], p),
        S.diff(ns["H"], x),
        S.diff(ns["H"], p)
    )
    for x, p in pairs
]

assert len(derivatives) == 11

print("PASS: 11 canonical pairs")
print("STAGE 1 PASS", flush=True)

def witness(w, rho, t):
    w = S.Rational(w)
    rho = S.Rational(rho)
    t = S.Rational(t)

    assert w != 0

    point = dict(ns["point"])

    point[q[1]] = w
    point[P[2]] = rho
    point[p1[3]] = t

    W0 = w**2
    z0 = rho**2

    IS0 = 12/W0 + S.Rational(448, 961)*z0

    point[ns["IS"]] = IS0

    IA0 = -S.Rational(3, 2)
    C0 = S.Rational(31, 8)

    M0 = 1 + (IA0**2 + IS0**2 + 1)/2

    point[p0[0]] = W0*z0/(4*C0) - M0

    return point

def evaluate(expr, point):
    return S.cancel(expr.subs(point))

def baseline(w, rho):
    W = w**2
    z = rho**2

    numerator = (
        476656*W**3*z
        + 802816*W**2*z**2
        + 19393941*W**2
        + 41330688*W*z
        + 531948096
    )

    return S.cancel(numerator/(14776336*W))

print("STAGE 2: exact transverse witnesses", flush=True)

cases = [
    (1, 1, 0),
    (1, 1, 1),
    (1, 1, 2),
    (2, 1, 1),
    (S.Rational(1, 2), 2, 1),
    (1, 0, 1)
]

results = {}

for w, rho, t in cases:
    w, rho, t = map(S.Rational, (w, rho, t))
    point = witness(w, rho, t)

    for name, expr in constraints:
        residual = evaluate(expr, point)
        assert residual == 0, (
            w, rho, t, name, residual
        )

    contributions = []

    for Bx, Bp, Hx, Hp in derivatives:
        term = S.cancel(
            evaluate(Bx, point)*evaluate(Hp, point)
            -
            evaluate(Bp, point)*evaluate(Hx, point)
        )
        contributions.append(term)

    D = S.cancel(sum(contributions))

    results[(w, rho, t)] = D

    print("w =", w, "rho =", rho, "t =", t)
    print("PASS: ten constraints")
    print("D =", D)
    print("sign =", S.sign(D))
    print("D minus t=0 baseline =",
          S.cancel(D-baseline(w, rho)),
          flush=True)

assert results[(S.Integer(1), S.Integer(1),
                S.Integer(0))] == baseline(
                    S.Integer(1), S.Integer(1)
                )

assert results[(S.Integer(1), S.Integer(1),
                S.Integer(1))] == (
                    baseline(S.Integer(1), S.Integer(1))
                    + S.Rational(8, 31)
                )

assert results[(S.Integer(1), S.Integer(1),
                S.Integer(2))] == (
                    baseline(S.Integer(1), S.Integer(1))
                    + S.Rational(32, 31)
                )

print("PASS: R28C and R27U regressions")
print("STAGE 2 PASS", flush=True)

print("STAGE 3: classification", flush=True)

for key, D in results.items():
    if D == 0:
        print("POTENTIAL D=0 WITNESS:", key)
    else:
        print("D != 0:", key)

print("STAGE 3 PASS", flush=True)

print()
print("R28D COMPLETE")
print("TRANSVERSE CANONICAL WITNESSES: CHECKED")
print("SYMBOLIC t DEPENDENCE: OPEN")
print("GENERAL IA BRANCH: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3)
)
