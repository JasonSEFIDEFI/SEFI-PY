import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29AN ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29AM.py").read_text()
marker = '\nprint("STAGE 1: reconstruct general cubic"'

assert source.count(marker) == 1, "STOP: source boundary ambiguous"

ns = {
    "__name__": "r29an_definitions",
    "__file__": str(folder / "CP098_T90C_R29AM.py")
}
exec(
    compile(source[:source.index(marker)], "R29AM definitions", "exec"),
    ns
)

x, y, h, b = ns["x"], ns["y"], ns["h"], ns["b"]
ND, NT = ns["ND"], ns["NT"]
reduce_expr = ns["reduce_expr"]

q = S.Symbol("q")

print("STAGE 1: reconstruct exact polynomials", flush=True)

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
Gq = S.Poly(G.as_expr(), q)
a0, a1, a2, a3 = [Gq.nth(k) for k in range(4)]

Delta = S.Poly(
    a1**2*a2**2
    - 4*a3*a1**3
    - 4*a0*a2**3
    - 27*a3**2*a0**2
    + 18*a3*a2*a1*a0,
    x, y
)

assert len(Delta.terms()) == 354

minimum_x = min(ex for (ex, ey), c in Delta.terms())
assert minimum_x == 4

D = S.Poly.from_dict(
    {(ex - 4, ey): c for (ex, ey), c in Delta.terms()},
    (x, y)
)

assert D.degree(x) == 26

print("PASS: exact B1 and reduced discriminant", flush=True)

print("STAGE 2: exact slice root counts", flush=True)

for uv in (0, 1, 4, 16):
    t0 = perf_counter()
    yv = S.Rational(31, 32) + uv

    B = S.Poly(a1.subs(y, yv), x)
    E = S.Poly(D.as_expr().subs(y, yv), x)

    assert B.degree() == 8
    assert E.degree() == 26
    assert B.eval(0) > 0
    assert B.LC() < 0
    assert E.eval(0) > 0
    assert E.LC() < 0

    nb = B.count_roots(0, S.oo)
    ne = E.count_roots(0, S.oo)

    print(f"SLICE u={uv}, y={yv}", flush=True)
    print("  B1 positive roots =", nb, flush=True)
    print("  DISC positive roots =", ne, flush=True)

    assert nb == 1, "STOP: B1 root-count regression"

    if ne != 1:
        print("  STOP-GATE: discriminant not uniquely rooted")
        print("  ROOT ORDERING: NOT ESTABLISHED")
        continue

    # Locate the unique roots with exact rational
    # endpoint signs and Sturm counts.
    def isolate(poly):
        lo = S.Integer(0)
        hi = S.Integer(1)

        assert poly.eval(lo) > 0

        while poly.eval(hi) >= 0:
            hi *= 2
            assert hi <= 2**20, (
                "STOP: root bracketing limit reached"
            )

        assert poly.count_roots(lo, hi) == 1

        for _ in range(24):
            mid = (lo + hi) / 2
            value = poly.eval(mid)

            if value == 0:
                return mid, mid

            if value > 0:
                lo = mid
            else:
                hi = mid

        assert poly.count_roots(lo, hi) == 1
        return lo, hi

    blo, bhi = isolate(B)
    dlo, dhi = isolate(E)

    print("  B1 root interval =", (blo, bhi), flush=True)
    print("  DISC root interval =", (dlo, dhi), flush=True)

    if dhi < blo:
        print("  PASS: DISC root precedes B1 root")
    elif bhi < dlo:
        print("  STOP-GATE: B1 root precedes DISC root")
    else:
        print("  STOP-GATE: root intervals overlap")

    print(
        "  elapsed seconds =",
        round(perf_counter() - t0, 3),
        flush=True
    )

print("R29AN COMPLETE")
print("RATIONAL-SLICE ROOT ORDERING: CHECKED")
print("CONTINUOUS-u ROOT ORDERING: OPEN")
print("GENERAL-y POSITIVITY: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter() - start, 3),
    flush=True
)
