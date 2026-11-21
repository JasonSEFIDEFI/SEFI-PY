import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29AJ ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29AI.py").read_text()

marker = '\nprint("STAGE 1: reconstruct exact shifted polynomial"'
assert source.count(marker) == 1, (
    "STOP: ambiguous R29AI source boundary"
)

ns = {
    "__name__": "r29aj_definitions",
    "__file__": str(folder / "CP098_T90C_R29AI.py")
}

exec(
    compile(source[:source.index(marker)],
            "R29AI definitions", "exec"),
    ns
)

x, y, h, b = ns["x"], ns["y"], ns["h"], ns["b"]
ND, NT = ns["ND"], ns["NT"]
reduce_expr = ns["reduce_expr"]

q = S.Symbol("q")

print("STAGE 1: reconstruct general cubic", flush=True)

PD = S.Poly(S.cancel(y * reduce_expr(ND)), b)
PT = S.Poly(S.cancel(y * reduce_expr(NT)), b)

R = S.Poly(
    S.resultant(PD.as_expr(), PT.as_expr(), b),
    x, y, h
)

assert len(R.terms()) == 63

terms = {}
for (ex, ey, eh), coeff in R.terms():
    assert ex >= 4
    assert eh % 2 == 0
    terms[(ex - 4, ey, eh // 2)] = coeff

G = S.Poly.from_dict(terms, (x, y, q))

assert G.degree(q) == 3

print("PASS: exact cubic reconstructed", flush=True)

print("STAGE 2: exact rational sample classification", flush=True)

positive_q_cases = []

for uv in (0, 1, 4, 16):
    yv = S.Rational(31, 32) + uv

    for xv in (8, 16, 32, 64):
        P = S.Poly(
            G.as_expr().subs({x: xv, y: yv}),
            q
        )

        assert P.degree() == 3

        a0, a1, a2, a3 = [
            P.nth(k) for k in range(4)
        ]

        assert a0 > 0
        assert a2 > 0
        assert a3 > 0

        delta = S.discriminant(P.as_expr(), q)

        b1sign = S.sign(a1)
        dsign = S.sign(delta)
        roots = P.count_roots(0, S.oo)

        assert b1sign in (-1, 0, 1)
        assert dsign in (-1, 0, 1)
        assert roots in (0, 1, 2)

        # Necessary and sufficient cubic criterion,
        # given a0,a2,a3 strictly positive.
        if a1 >= 0:
            assert roots == 0
        elif delta < 0:
            assert roots == 0
        elif delta > 0:
            assert roots == 2
        else:
            assert roots == 1

        print(
            f"u={uv:2d} x={xv:2d} "
            f"B1={str(b1sign):>2s} "
            f"DISC={str(dsign):>2s} "
            f"POSITIVE_Q_ROOTS={roots}",
            flush=True
        )

        if roots:
            positive_q_cases.append((uv, xv, roots))

print("STAGE 3: classification", flush=True)

print(
    "POSITIVE-q SAMPLE CASES =",
    positive_q_cases,
    flush=True
)

if positive_q_cases:
    print(
        "STOP-GATE: positive resultant roots detected",
        flush=True
    )
    print(
        "ORIGINAL EQUATION RECOVERY REQUIRED",
        flush=True
    )
else:
    print(
        "PASS: no positive resultant roots "
        "at tested rational points",
        flush=True
    )

print("R29AJ COMPLETE")
print("CONTINUOUS GENERAL-y DOMAIN: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter() - start, 3),
    flush=True
)
