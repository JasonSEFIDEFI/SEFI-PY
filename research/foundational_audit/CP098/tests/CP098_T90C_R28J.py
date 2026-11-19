import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28J ===", flush=True)

source = Path(__file__).with_name("CP098_T90C_R27U.py")
assert source.is_file()

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
assert len(pairs) == 11

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

print("STAGE 1: independent canonical derivatives",
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

print("PASS: 11 canonical pairs", flush=True)
print("STAGE 1 PASS", flush=True)

base = S.Rational(593952197, 14776336)

def witness(z, b):
    point = dict(ns["point"])

    point[q[1]] = S.Integer(1)
    point[P[2]] = S.Integer(1)
    point[p1[2]] = z
    point[p1[3]] = S.Integer(0)

    point[ns["pA"]] = S.Integer(0)
    point[ns["pS"]] = S.Integer(0)
    point[ns["pW"]] = b

    IS = S.Integer(12) + S.Rational(448, 961)
    point[ns["IS"]] = IS

    M = 1 + (
        S.Rational(9, 4) + IS**2 + 1 + b**2
    )/2

    point[p0[0]] = S.Rational(2, 31) - M
    point[p0[1]] = -S.Rational(8, 31)*z

    return point

def at(expr, point):
    return S.simplify(expr.subs(point))

def check(label, z, b, expected):
    print("WITNESS:", label, flush=True)

    point = witness(z, b)

    for name, expr in constraints:
        residual = at(expr, point)
        assert residual == 0, (
            label, name, residual
        )

    print("PASS: ten constraints", flush=True)

    contributions = []
    for Bx, Bp, Hx, Hp in derivatives:
        term = S.simplify(
            at(Bx, point)*at(Hp, point)
            - at(Bp, point)*at(Hx, point)
        )
        contributions.append(term)

    D = S.simplify(sum(contributions))

    predicted = S.simplify(
        base
        + S.Rational(8, 31)*z**2
        + b**2/4
        + S.Rational(16, 31)*b*z
    )

    assert S.simplify(D-predicted) == 0
    assert S.simplify(D-expected) == 0

    print("PASS: independent bracket regression")
    print("D =", D, flush=True)

    return D

print("STAGE 2: sign-changing rational witnesses",
      flush=True)

positive = check(
    "positive baseline",
    S.Integer(0),
    S.Integer(0),
    base
)

z_negative = S.Integer(100)
b_negative = -S.Rational(32, 31)*z_negative

negative_expected = S.simplify(
    base - S.Rational(8, 961)*z_negative**2
)

negative = check(
    "negative rational witness",
    z_negative,
    b_negative,
    negative_expected
)

assert positive > 0
assert negative < 0

print("PASS: both signs realized", flush=True)
print("STAGE 2 PASS", flush=True)

print("STAGE 3: exact D=0 witness",
      flush=True)

z_star = S.sqrt(S.Rational(593952197, 123008))
b_star = -S.Rational(32, 31)*z_star

assert S.simplify(
    base - S.Rational(8, 961)*z_star**2
) == 0

zero = check(
    "exact algebraic zero",
    z_star,
    b_star,
    S.Integer(0)
)

assert zero == 0

print("z_star =", z_star, flush=True)
print("b_star =", b_star, flush=True)
print("PASS: exact canonical D=0", flush=True)
print("STAGE 3 PASS", flush=True)

print()
print("R28J COMPLETE")
print("TEN CONSTRAINTS: VERIFIED AT ALL WITNESSES")
print("POSITIVE AND NEGATIVE D: VERIFIED")
print("EXACT D=0 WITNESS: VERIFIED")
print("NEXT DIRAC PRESERVATION: NOT YET CHECKED")
print("GENERAL IA BRANCH: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3)
)
