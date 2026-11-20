import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28F ===", flush=True)

source = Path(__file__).with_name(
    "CP098_T90C_R27U.py"
)
assert source.is_file(), "R27U source missing"

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

rho, z = S.symbols("rho z", real=True)

u, q, P = ns["u"], ns["q"], ns["P"]
p0, p1 = ns["p0"], ns["p1"]

point = dict(ns["point"])

point[P[2]] = rho
point[p1[2]] = z
point[p1[3]] = S.Integer(0)

IS0 = 12 + S.Rational(448, 961)*rho**2
point[ns["IS"]] = IS0

IA0 = -S.Rational(3, 2)
C0 = S.Rational(31, 8)

M0 = 1 + (IA0**2 + IS0**2 + 1)/2

# Preserve theta=0.
point[p0[0]] = rho**2/(4*C0) - M0

# Preserve B=0 when p1 is parallel to P.
point[p0[1]] = -S.Rational(8, 31)*rho*z

def at(expr):
    return S.cancel(expr.subs(point))

print("STAGE 1: symbolic constraint audit",
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

print("STAGE 2: canonical D(rho,z)",
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

assert len(pairs) == 11

parts = []

for i, (x, p) in enumerate(pairs, 1):
    bx = at(S.diff(ns["B"], x))
    bp = at(S.diff(ns["B"], p))
    hx = at(S.diff(ns["H"], x))
    hp = at(S.diff(ns["H"], p))

    term = S.cancel(bx*hp - bp*hx)
    parts.append(term)

    print(f"pair {i:02d} = {term}",
          flush=True)

D = S.factor(sum(parts))

print("D(rho,z) =", D, flush=True)
print("STAGE 2 PASS", flush=True)

print("STAGE 3: independent regressions",
      flush=True)

R27Y = (
    S.Rational(50176, 923521)*rho**4
    + S.Rational(2719, 961)*rho**2
    + S.Rational(597, 16)
)

assert S.cancel(D.subs(z, 0) - R27Y) == 0
print("PASS: R27Y recovered at z=0")

baseline = S.Rational(
    593952197, 14776336
)

assert S.cancel(
    D.subs({rho: 1, z: 0}) - baseline
) == 0

# R27V established this exact z correction
# on the rho=1, W=1 family.
expected_r27v = baseline + S.Rational(8, 31)*z**2

assert S.cancel(
    D.subs(rho, 1) - expected_r27v
) == 0

print("PASS: R27V z-family recovered")
print("STAGE 3 PASS", flush=True)

print("STAGE 4: new coupling analysis",
      flush=True)

correction = S.factor(D - R27Y)

print("D minus R27Y =", correction)

# Classify the z dependence exactly.
poly_z = S.Poly(D, z)

print("degree in z =", poly_z.degree())
print("coefficients =", poly_z.all_coeffs())

if poly_z.degree() == 2:
    a = poly_z.coeff_monomial(z**2)
    b = poly_z.coeff_monomial(z)
    c = poly_z.coeff_monomial(1)

    discriminant = S.factor(b**2 - 4*a*c)

    print("quadratic a =", a)
    print("quadratic b =", b)
    print("quadratic c =", c)
    print("discriminant =", discriminant)

    if S.simplify(a > 0) is S.true:
        print("leading coefficient: positive")
    else:
        print("leading coefficient: requires analysis")

    print(
        "NOTE: discriminant and sign require "
        "exact interpretation over real rho"
    )
else:
    print(
        "Quadratic classification not applicable"
    )

print("STAGE 4 COMPLETE", flush=True)

print()
print("R28F COMPLETE")
print("PARALLEL MOMENTUM FAMILY: DERIVED")
print("R27Y/R27V REGRESSIONS: CHECKED")
print("GLOBAL SIGN: REQUIRES INTERPRETATION")
print("GENERAL IA BRANCH: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter() - start, 3)
)
