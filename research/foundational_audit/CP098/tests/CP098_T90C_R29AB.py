import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29AB ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29AA.py").read_text()
cut = source.index('print("STAGE 1: reconstruct exact certificate"')

ns = {
    "__name__": "r29ab_definitions",
    "__file__": str(folder / "CP098_T90C_R29AA.py")
}
exec(compile(source[:cut], "R29AA definitions", "exec"), ns)

x, y, h, b, q = (
    ns["x"], ns["y"], ns["h"], ns["b"], ns["q"]
)
ND, NT = ns["ND"], ns["NT"]
reduce_expr = ns["reduce_expr"]

print("STAGE 1: reconstruct cubic", flush=True)

PD = S.cancel(y * reduce_expr(ND))
PT = S.cancel(y * reduce_expr(NT))

A = S.expand(PD.subs(y, 4))
B = S.expand(PT.subs(y, 4))

R = S.Poly(S.resultant(A, B, b), x, h)

assert len(R.terms()) == 24
assert R.degree(x) == 13
assert R.degree(h) == 6

terms = {}
for (ex, eh), coeff in R.terms():
    assert ex >= 4 and eh % 2 == 0
    terms[(ex - 4, eh // 2)] = coeff

F = S.Poly.from_dict(terms, (x, q))
assert F.degree(q) == 3

Fq = S.Poly(F.as_expr(), q)
a0, a1, a2, a3 = [Fq.nth(i) for i in range(4)]

for coeff in (a0, a2, a3):
    assert all(
        c.is_positive is True
        for c in S.Poly(coeff, x).coeffs()
    )

print("PASS: exact cubic and positive A0,A2,A3", flush=True)

print("STAGE 2: cubic discriminant", flush=True)

# Explicit cubic discriminant avoids a generic
# symbolic resultant or expensive factorization.
disc = S.Poly(
    a1**2 * a2**2
    - 4*a3*a1**3
    - 4*a0*a2**3
    - 27*a3**2*a0**2
    + 18*a3*a2*a1*a0,
    x
)

assert disc.as_expr() != 0

print("DISCRIMINANT DEGREE =", disc.degree(), flush=True)
print("DISCRIMINANT TERMS =", len(disc.terms()), flush=True)

positive = sum(
    1 for c in disc.coeffs()
    if c.is_positive is True
)
negative = sum(
    1 for c in disc.coeffs()
    if c.is_negative is True
)

assert positive + negative == len(disc.coeffs())

print("POSITIVE COEFFICIENTS =", positive, flush=True)
print("NEGATIVE COEFFICIENTS =", negative, flush=True)

print("STAGE 3: exact rational sign checks", flush=True)

values = [
    S.Rational(1, 100),
    S.Rational(1, 10),
    S.Integer(1),
    S.Integer(10),
    S.Integer(100)
]

for value in values:
    result = disc.eval(value)
    assert result.is_Rational
    sign = S.sign(result)
    assert sign in (S.Integer(-1), S.Integer(0), S.Integer(1))
    print("x =", value, "DISC SIGN =", sign, flush=True)

print("STAGE 4: certificate assessment", flush=True)

if positive == 0:
    print("CERTIFIED: discriminant < 0 for x > 0")
    print("CERTIFIED: F(x,q) > 0 for x > 0, q > 0")
    print("CERTIFIED: no simultaneous zeros on y=4")
else:
    print("GLOBAL DISCRIMINANT SIGN: OPEN")
    print("FULL CUBIC POSITIVITY: OPEN")
    print("No conclusion from sampled signs alone")

print("R29AB COMPLETE")
print("GENERAL y DOMAIN: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3),
    flush=True
)
