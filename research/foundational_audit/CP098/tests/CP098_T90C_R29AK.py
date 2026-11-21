import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29AK ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29AJ.py").read_text()
marker = '\nprint("STAGE 1: reconstruct general cubic"'

assert source.count(marker) == 1, (
    "STOP: ambiguous source boundary"
)

ns = {
    "__name__": "r29ak_definitions",
    "__file__": str(folder / "CP098_T90C_R29AJ.py")
}

exec(
    compile(source[:source.index(marker)],
            "R29AJ definitions", "exec"),
    ns
)

x, y, h, b = (
    ns["x"], ns["y"], ns["h"], ns["b"]
)
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

assert len(B1.terms()) == 84
assert B1.degree(x) == 8

print("PASS: B1 reconstructed", flush=True)

print("STAGE 2: exact coefficient signs", flush=True)

Px = S.Poly(B1.as_expr(), x)
fixed_signs = []

for k in range(9):
    C = S.Poly(Px.nth(k), u)

    positives = sum(1 for c in C.coeffs() if c.is_positive is True)
    negatives = sum(1 for c in C.coeffs() if c.is_negative is True)

    assert positives + negatives == len(C.coeffs()), (
        f"STOP: unresolved coefficient sign at x^{k}"
    )

    if positives and not negatives:
        sign = 1
    elif negatives and not positives:
        sign = -1
    elif not positives and not negatives:
        sign = 0
    else:
        sign = None

    fixed_signs.append(sign)

    print(
        f"x^{k}: degree_u={C.degree()} "
        f"positive={positives} "
        f"negative={negatives} "
        f"fixed_sign={sign}",
        flush=True
    )

print("STAGE 3: Descartes certificate", flush=True)

if any(s is None for s in fixed_signs):
    print("STOP-GATE: mixed-sign coefficient polynomial")
    print("UNIVERSAL ROOT COUNT: OPEN")

else:
    signs = [s for s in fixed_signs if s != 0]

    changes = sum(
        signs[i] != signs[i - 1]
        for i in range(1, len(signs))
    )

    print("SIGN CHANGES =", changes, flush=True)

    if signs[0] == 1 and signs[-1] == -1 and changes == 1:
        print(
            "CERTIFIED: exactly one positive "
            "B1 root for every u>0",
            flush=True
        )
        print(
            "CERTIFIED: B1 positive before "
            "root and negative after",
            flush=True
        )
    else:
        print("STOP-GATE: direct Descartes certificate fails")
        print("UNIVERSAL ROOT COUNT: OPEN")

print("R29AK COMPLETE")
print("GENERAL-y CUBIC POSITIVITY: OPEN")
print("DISCRIMINANT SIGN OVERLAP: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter() - start, 3),
    flush=True
)
