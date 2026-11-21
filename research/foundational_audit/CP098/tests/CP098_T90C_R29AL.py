import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29AL ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29AK.py").read_text()

marker = '\nprint("STAGE 1: reconstruct B1"'
assert source.count(marker) == 1, "STOP: source boundary ambiguous"

ns = {
    "__name__": "r29al_definitions",
    "__file__": str(folder / "CP098_T90C_R29AK.py")
}

exec(
    compile(source[:source.index(marker)],
            "R29AK definitions", "exec"),
    ns
)

x, y, h, b = ns["x"], ns["y"], ns["h"], ns["b"]
ND, NT = ns["ND"], ns["NT"]
reduce_expr = ns["reduce_expr"]

u, q = S.symbols("u q")

print("STAGE 1: reconstruct B1", flush=True)

PD = S.Poly(S.cancel(y * reduce_expr(ND)), b)
PT = S.Poly(S.cancel(y * reduce_expr(NT)), b)

R = S.Poly(
    S.resultant(PD.as_expr(), PT.as_expr(), b),
    x, y, h
)

assert len(R.terms()) == 63

terms = {}
for (ex, ey, eh), coeff in R.terms():
    assert ex >= 4 and eh % 2 == 0
    terms[(ex - 4, ey, eh // 2)] = coeff

G = S.Poly.from_dict(terms, (x, y, q))

H = S.Poly(
    G.as_expr().subs(y, u + S.Rational(31, 32)),
    x, u, q
)

assert len(H.terms()) == 203

B1 = S.Poly(
    S.Poly(H.as_expr(), q).nth(1),
    x, u
)

assert B1.degree(x) == 8
assert len(B1.terms()) == 84

print("PASS: B1 reconstructed", flush=True)

print("STAGE 2: verify fixed-sign coefficients", flush=True)

Px = S.Poly(B1.as_expr(), x)

for k in range(9):
    C = S.Poly(Px.nth(k), u)
    coeffs = C.coeffs()

    if k <= 5:
        assert coeffs
        assert all(c.is_positive is True for c in coeffs)
        print(f"PASS: x^{k} coefficient positive for u>0")

    elif k == 6:
        assert coeffs
        assert all(
            c.is_positive is True or c.is_negative is True
            for c in coeffs
        )
        print("PASS: x^6 coefficient unrestricted")

    else:
        assert coeffs
        assert all(c.is_negative is True for c in coeffs)
        print(f"PASS: x^{k} coefficient negative for u>0")

print("STAGE 3: universal Descartes argument", flush=True)

# Descending x coefficients:
# x^8,x^7: negative
# x^6: arbitrary sign or zero
# x^5,...,x^0: positive
#
# Every possible sign sequence has exactly one
# sign change after zero coefficients are removed.

for middle in (-1, 0, 1):
    signs = [-1, -1, middle] + [1] * 6
    nonzero = [s for s in signs if s != 0]

    changes = sum(
        nonzero[i] != nonzero[i - 1]
        for i in range(1, len(nonzero))
    )

    assert changes == 1

print("PASS: exactly one Descartes sign change")
print("PASS: at most one positive root")

# For every u>0, B1(0,u)>0 and its
# leading x coefficient is negative.
# Therefore B1(x,u) tends to -infinity
# as x tends to +infinity.

print("PASS: opposite endpoint signs")
print("PASS: at least one positive root")

print("STAGE 4: exact conclusion", flush=True)

print("CERTIFIED: exactly one positive B1 root for every u>0")
print("CERTIFIED: B1 positive before that root")
print("CERTIFIED: B1 negative after that root")

print("R29AL COMPLETE")
print("UNIVERSAL B1 ROOT COUNT: VERIFIED")
print("DISCRIMINANT SIGN OVERLAP: OPEN")
print("GENERAL-y CUBIC POSITIVITY: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter() - start, 3),
    flush=True
)
