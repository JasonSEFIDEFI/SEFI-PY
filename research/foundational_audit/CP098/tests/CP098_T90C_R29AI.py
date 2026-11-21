import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29AI ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29AH.py").read_text()
marker = '\nprint("STAGE 1: reconstruct general resultant"'

assert source.count(marker) == 1, "STOP: source boundary ambiguous"

ns = {
    "__name__": "r29ai_definitions",
    "__file__": str(folder / "CP098_T90C_R29AH.py")
}
exec(compile(source[:source.index(marker)], "R29AH definitions", "exec"), ns)

x, y, h, b = ns["x"], ns["y"], ns["h"], ns["b"]
ND, NT = ns["ND"], ns["NT"]
reduce_expr = ns["reduce_expr"]

u = S.Symbol("u")
q = S.Symbol("q")

print("STAGE 1: reconstruct exact shifted polynomial", flush=True)

PD = S.Poly(S.cancel(y * reduce_expr(ND)), b)
PT = S.Poly(S.cancel(y * reduce_expr(NT)), b)

R = S.Poly(S.resultant(PD.as_expr(), PT.as_expr(), b), x, y, h)
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

B1 = S.Poly(S.Poly(H.as_expr(), q).nth(1), x, u)
assert len(B1.terms()) == 84

print("PASS: R29AH B1 reconstructed", flush=True)

print("STAGE 2: exact positive-root counts", flush=True)

for uv in [0, 1, 4, 16]:
    P = S.Poly(B1.as_expr().subs(u, uv), x)
    assert not P.is_zero

    count = P.count_roots(0, S.oo)

    print("u =", uv, "y =", S.Rational(31, 32) + uv)
    print("  degree =", P.degree())
    print("  positive roots =", count)
    print("  sign at x=0 =", S.sign(P.eval(0)))
    print("  leading sign =", S.sign(P.LC()))

    for xv in [S.Rational(1, 10), 1, 4, 8, 16, 32]:
        sign = S.sign(P.eval(xv))
        assert sign in (-1, 0, 1)
        print("  x =", xv, "B1 sign =", sign)

    print("  exact slice complete", flush=True)

print("R29AI COMPLETE")
print("B1 SIGN STRUCTURE: EXACT SLICES ONLY")
print("GENERAL-y POSITIVITY: OPEN")
print("GENERAL-y SIMULTANEOUS ZERO CLASSIFICATION: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =", round(perf_counter()-start, 3))
