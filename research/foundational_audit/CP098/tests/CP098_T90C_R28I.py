import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28I ===", flush=True)

source = Path(__file__).with_name(
    "CP098_T90C_R27U.py"
)
assert source.is_file(), "R27U source missing"

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

w = S.symbols("w", real=True, nonzero=True)
rho, z, t, a, b = S.symbols(
    "rho z t a b", real=True
)

u, q, P = ns["u"], ns["q"], ns["P"]
p0, p1 = ns["p0"], ns["p1"]

point = dict(ns["point"])

point[q[1]] = w
point[P[2]] = rho
point[p1[2]] = z
point[p1[3]] = t
point[ns["pS"]] = a
point[ns["pW"]] = b

W0 = w**2
r2 = rho**2
IS0 = 12/W0 + S.Rational(448, 961)*r2
point[ns["IS"]] = IS0

M0 = (
    1
    + (
        S.Rational(9, 4)
        + IS0**2
        + 1
        + a**2
        + b**2
    )/2
)

point[p0[0]] = W0*r2/S.Rational(31, 2) - M0
point[p0[1]] = -S.Rational(8, 31)*w*rho*z

def at(expr):
    return S.cancel(expr.subs(point))

print("STAGE 1: symbolic constraints", flush=True)

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
    print("PASS:", name, flush=True)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: full canonical bracket", flush=True)

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

parts = []

for i, (x, p) in enumerate(pairs, 1):
    # Differentiate BEFORE imposing the constraints.
    Bx = at(S.diff(ns["B"], x))
    Bp = at(S.diff(ns["B"], p))
    Hx = at(S.diff(ns["H"], x))
    Hp = at(S.diff(ns["H"], p))

    term = S.cancel(Bx*Hp - Bp*Hx)
    parts.append(term)

    print(f"pair {i:02d} = {term}", flush=True)

D = S.cancel(sum(parts))

print("D derived", flush=True)
print("STAGE 2 PASS", flush=True)

print("STAGE 3: isolate identity correction", flush=True)

baseline = (
    476656*W0**3*r2
    + 802816*W0**2*r2**2
    + 19393941*W0**2
    + 41330688*W0*r2
    + 531948096
)/(14776336*W0) + (
    S.Rational(8, 31)*W0*(z**2+t**2)
)

assert S.cancel(
    D.subs({a: 0, b: 0}) - baseline
) == 0

print("PASS: R28G recovered", flush=True)

correction = S.factor(D - baseline)
print("identity correction =", correction, flush=True)

poly = S.Poly(correction, a, b)

print("identity-momentum total degree =",
      poly.total_degree(), flush=True)
print("monomials and coefficients:", flush=True)

for powers, coeff in poly.terms():
    print(powers, ":", S.factor(coeff), flush=True)

print("STAGE 3 PASS", flush=True)

print("STAGE 4: exact witness regressions", flush=True)

checks = [
    (1, 1, 0, 0, 1, 0, S.Rational(1, 4)),
    (1, 1, 0, 0, 0, 2, S.Integer(1)),
    (1, 1, 0, 0, 1, 2, S.Rational(5, 4)),
    (2, 1, 1, 1, 1, 1, S.Rational(126, 31)),
    (
        S.Rational(1, 2), 2, 1, 0, 1, 2,
        S.Rational(411, 496)
    )
]

for wv, rv, zv, tv, av, bv, expected in checks:
    actual = S.cancel(correction.subs({
        w: wv, rho: rv, z: zv, t: tv,
        a: av, b: bv
    }))
    assert actual == expected, (
        (wv, rv, zv, tv, av, bv),
        actual, expected
    )
    print("PASS:", (wv, rv, zv, tv, av, bv),
          flush=True)

print("STAGE 4 PASS", flush=True)

print("STAGE 5: sign diagnostics", flush=True)

if poly.total_degree() == 2:
    aa = S.factor(poly.coeff_monomial(a**2))
    ab = S.factor(poly.coeff_monomial(a*b))
    bb = S.factor(poly.coeff_monomial(b**2))
    la = S.factor(poly.coeff_monomial(a))
    lb = S.factor(poly.coeff_monomial(b))
    c = S.factor(poly.coeff_monomial(1))

    print("a² coefficient =", aa)
    print("ab coefficient =", ab)
    print("b² coefficient =", bb)
    print("a coefficient =", la)
    print("b coefficient =", lb)
    print("constant =", c)
    print("quadratic determinant =",
          S.factor(4*aa*bb-ab**2))
else:
    print("Nonquadratic correction: investigate")

print("STAGE 5 COMPLETE", flush=True)

print()
print("R28I COMPLETE")
print("SIX-VARIABLE CANONICAL FAMILY: DERIVED")
print("R28G/R28H REGRESSIONS: VERIFIED")
print("GLOBAL SIGN: PENDING EXACT INTERPRETATION")
print("GENERAL IA BRANCH: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3)
)
