import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29AD ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29AC.py").read_text()
cut = source.index('print("STAGE 1: reconstruct exact cubic"')

ns = {
    "__name__": "r29ad_definitions",
    "__file__": str(folder / "CP098_T90C_R29AC.py")
}
exec(compile(source[:cut], "R29AC definitions", "exec"), ns)

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
for (ex, eh), coefficient in R.terms():
    assert ex >= 4 and eh % 2 == 0
    terms[(ex - 4, eh // 2)] = coefficient

F = S.Poly.from_dict(terms, (x, q))
Fq = S.Poly(F.as_expr(), q)

a0, a1, a2, a3 = [Fq.nth(i) for i in range(4)]

A1 = S.Poly(a1, x)

disc = S.Poly(
    a1**2*a2**2
    - 4*a3*a1**3
    - 4*a0*a2**3
    - 27*a3**2*a0**2
    + 18*a3*a2*a1*a0,
    x
)

assert A1.degree() == 8
assert disc.degree() == 30

# Remove only the exact x=0 boundary factor.
# This is valid because the admissible domain is x>0.
disc_terms = disc.terms()
minimum_power = min(
    powers[0] for powers, coefficient in disc_terms
)

assert minimum_power >= 1, (
    "STOP: expected boundary factor absent"
)

disc_reduced = S.Poly.from_dict(
    {
        (powers[0] - minimum_power,): coefficient
        for powers, coefficient in disc_terms
    },
    x
)

assert S.expand(
    disc.as_expr()
    - x**minimum_power * disc_reduced.as_expr()
) == 0

assert disc_reduced.eval(0) != 0

print(
    "DISCRIMINANT x=0 MULTIPLICITY =",
    minimum_power,
    flush=True
)

print(
    "REDUCED DISCRIMINANT DEGREE =",
    disc_reduced.degree(),
    flush=True
)

print("PASS: exact boundary-factor removal", flush=True)
print("PASS: R29AC polynomials reconstructed", flush=True)

print("STAGE 2: exact Sturm root counts", flush=True)

intervals = [
    ("(0,1)", S.Integer(0), S.Integer(1)),
    ("(1,10)", S.Integer(1), S.Integer(10)),
    ("(10,infinity)", S.Integer(10), S.oo),
    ("(0,infinity)", S.Integer(0), S.oo),
]

for name, poly in [
    ("A1", A1),
    ("DISC_REDUCED", disc_reduced)
]:
    print("POLYNOMIAL:", name, flush=True)

    counts = {}

    for label, left, right in intervals:
        # count_roots includes endpoint roots. Check
        # boundaries explicitly to preserve interpretation.
        if left != S.oo:
            assert poly.eval(left) != 0, (
                f"STOP: {name} root at {left}"
            )
        if right != S.oo:
            assert poly.eval(right) != 0, (
                f"STOP: {name} root at {right}"
            )

        count = poly.count_roots(left, right)
        counts[label] = count

        print(
            "  ROOTS", label, "=", count,
            flush=True
        )

    assert (
        counts["(0,1)"]
        + counts["(1,10)"]
        + counts["(10,infinity)"]
        == counts["(0,infinity)"]
    ), f"STOP: {name} interval count mismatch"

    print("PASS:", name, "root counts consistent", flush=True)

print("R29AD COMPLETE")
print("EXACT POSITIVE ROOT COUNTS: VERIFIED")
print("ROOT ORDERING AND SIGN OVERLAP: OPEN")
print("GLOBAL y=4 CLASSIFICATION: OPEN")
print("GENERAL y DOMAIN: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3),
    flush=True
)
