import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27Y ===", flush=True)

source = Path(__file__).with_name(
    "CP098_T90C_R27U.py"
)
assert source.is_file(), "R27U source missing"

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

rho = S.symbols("rho", real=True)
point = dict(ns["point"])

u, q, P = ns["u"], ns["q"], ns["P"]
p0, p1 = ns["p0"], ns["p1"]

point[p1[3]] = 0
point[P[2]] = rho

IS0 = 12 + S.Rational(448, 961)*rho**2
point[ns["IS"]] = IS0

IA0 = -S.Rational(3, 2)
C0 = S.Rational(31, 8)
M0 = 1 + (IA0**2 + IS0**2 + 1)/2

point[p0[0]] = rho**2/(4*C0) - M0

def at(expr):
    return S.cancel(expr.subs(point))

print("STAGE 1: symbolic constraint family",
      flush=True)

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

print("STAGE 2: canonical D(rho)",
      flush=True)

pairs = (
    list(zip(u, p1))
    + list(zip(q, P))
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

expected = {
    0: S.Rational(597, 16),
    1: S.Rational(593952197, 14776336),
    2: S.Rational(731416469, 14776336)
}

for value, target in expected.items():
    result = S.cancel(D.subs(rho, value))
    assert result == target, (value, result, target)
    print(f"PASS: rho={value} matches R27X")

print("D(rho) =", D)
print("STAGE 2 PASS", flush=True)

print("STAGE 3: polynomial sign analysis",
      flush=True)

num, den = S.fraction(S.cancel(D))
poly = S.Poly(num, rho)

print("numerator =", S.factor(num))
print("denominator =", den)
print("degree =", poly.degree())
print("coefficients =", poly.all_coeffs())

assert not poly.is_zero

roots = S.solveset(
    S.Eq(D, 0),
    rho,
    domain=S.S.Reals
)

print("real zeros =", roots)

# For a polynomial in rho^2, positive
# coefficients give an immediate certificate.
p = S.symbols("p", nonnegative=True)

if all(
    power[0] % 2 == 0
    for power, coeff in poly.terms()
):
    even_poly = S.Poly.from_dict(
        {
            (power[0]//2,): coeff
            for power, coeff in poly.terms()
        },
        p
    )

    print("numerator as polynomial in rho² =",
          even_poly.as_expr())

    if den.is_positive and all(
        coeff > 0
        for coeff in even_poly.all_coeffs()
    ):
        print(
            "PASS: D(rho)>0 for every real rho"
        )
    else:
        print(
            "SIGN CERTIFICATE: NOT ESTABLISHED"
        )
else:
    print(
        "SIGN CERTIFICATE: EVEN-POLYNOMIAL "
        "TEST NOT APPLICABLE"
    )

print("STAGE 3 PASS", flush=True)

print()
print("R27Y COMPLETE")
print("SYMBOLIC GEOMETRIC FAMILY: CHECKED")
print("REAL D=0 ROOTS: REPORTED")
print("GENERAL IA BRANCH: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3)
)
