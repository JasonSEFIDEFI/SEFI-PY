import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29AC ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29AB.py").read_text()
cut = source.index('print("STAGE 1: reconstruct cubic"')

ns = {
    "__name__": "r29ac_definitions",
    "__file__": str(folder / "CP098_T90C_R29AB.py")
}

exec(compile(source[:cut], "R29AB definitions", "exec"), ns)

x, y, h, b, q = (
    ns["x"], ns["y"], ns["h"], ns["b"], ns["q"]
)
ND, NT = ns["ND"], ns["NT"]
reduce_expr = ns["reduce_expr"]

print("STAGE 1: reconstruct exact cubic", flush=True)

PD = S.cancel(y * reduce_expr(ND))
PT = S.cancel(y * reduce_expr(NT))

A = S.expand(PD.subs(y, 4))
B = S.expand(PT.subs(y, 4))

R = S.Poly(S.resultant(A, B, b), x, h)

assert len(R.terms()) == 24

terms = {}
for (ex, eh), coeff in R.terms():
    assert ex >= 4 and eh % 2 == 0
    terms[(ex - 4, eh // 2)] = coeff

F = S.Poly.from_dict(terms, (x, q))
Fq = S.Poly(F.as_expr(), q)

a0, a1, a2, a3 = [Fq.nth(i) for i in range(4)]

for coeff in (a0, a2, a3):
    assert all(
        c.is_positive is True
        for c in S.Poly(coeff, x).coeffs()
    )

print("PASS: cubic structure verified", flush=True)

print("STAGE 2: reconstruct discriminant", flush=True)

disc = S.Poly(
    a1**2*a2**2
    - 4*a3*a1**3
    - 4*a0*a2**3
    - 27*a3**2*a0**2
    + 18*a3*a2*a1*a0,
    x
)

assert disc.degree() == 30

print("PASS: discriminant reconstructed", flush=True)

print("STAGE 3: exact positive-root tests", flush=True)

values = [
    S.Rational(1, 100),
    S.Rational(1, 10),
    S.Integer(1),
    S.Integer(10),
    S.Integer(100)
]

for value in values:
    a1v = S.Poly(a1, x).eval(value)
    dv = disc.eval(value)

    assert a1v.is_Rational
    assert dv.is_Rational

    a1_sign = S.sign(a1v)
    disc_sign = S.sign(dv)

    assert a1_sign in (-1, 0, 1)
    assert disc_sign in (-1, 0, 1)

    # Under a0,a2,a3 > 0:
    # a1 >= 0 excludes positive roots.
    # a1 < 0 and disc > 0 gives two positive roots.
    # a1 < 0 and disc = 0 gives a positive double root.
    expected = (
        2 if a1_sign == -1 and disc_sign == 1
        else 1 if a1_sign == -1 and disc_sign == 0
        else 0
    )

    specialized = S.Poly(
        F.as_expr().subs(x, value), q
    )

    actual = specialized.count_roots(0, S.oo)

    assert actual == expected, (
        f"STOP: root-count mismatch at x={value}: "
        f"predicted={expected}, actual={actual}"
    )

    print(
        "x =", value,
        "A1 sign =", a1_sign,
        "DISC sign =", disc_sign,
        "POSITIVE DISTINCT ROOTS =", actual,
        flush=True
    )

print("PASS: exact root classifications agree", flush=True)

print("R29AC COMPLETE")
print("RATIONAL SLICES: CLASSIFIED")
print("GLOBAL y=4 CLASSIFICATION: OPEN")
print("GENERAL y DOMAIN: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter() - start, 3),
    flush=True
)
