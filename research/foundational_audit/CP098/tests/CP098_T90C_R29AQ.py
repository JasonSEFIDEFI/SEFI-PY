import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29AQ ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29AP.py").read_text()

marker = '\nprint("STAGE 1: reconstruct exact discriminant"'
assert source.count(marker) == 1, "STOP: source boundary ambiguous"

ns = {
    "__name__": "r29aq_definitions",
    "__file__": str(folder / "CP098_T90C_R29AP.py")
}

exec(
    compile(source[:source.index(marker)],
            "R29AP definitions", "exec"),
    ns
)

x, y, h, b = ns["x"], ns["y"], ns["h"], ns["b"]
ND, NT = ns["ND"], ns["NT"]
reduce_expr = ns["reduce_expr"]

q, u = S.symbols("q u")

print("STAGE 1: reconstruct discriminant", flush=True)

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

min_x = min(ex for (ex, ey), c in Delta.terms())
assert min_x == 4

D = S.Poly.from_dict(
    {(ex-min_x, ey): c
     for (ex, ey), c in Delta.terms()},
    (x, y)
)

shifted = S.Poly(
    D.as_expr().subs(y, u + S.Rational(31,32)),
    x, u
)
assert len(shifted.terms()) == 833

Px = S.Poly(shifted.as_expr(), x)

print("PASS: discriminant reconstructed", flush=True)

print("STAGE 2: isolate coefficient roots", flush=True)

intervals = {}

for k in range(14, 23):
    C = S.Poly(Px.nth(k), u)

    assert C.eval(0) < 0
    assert C.LC() > 0

    # Recheck the exact Descartes certificate.
    signs = [
        1 if c > 0 else -1 if c < 0 else 0
        for c in C.all_coeffs()
    ]
    signs = [s for s in signs if s]

    changes = sum(
        signs[i] != signs[i-1]
        for i in range(1, len(signs))
    )
    assert changes == 1

    lo = S.Integer(0)
    hi = S.Integer(1)

    while C.eval(hi) <= 0:
        hi *= 2
        assert hi <= 2**20, (
            "STOP: root bracketing limit"
        )

    for _ in range(48):
        mid = (lo + hi) / 2
        value = C.eval(mid)

        if value == 0:
            lo = hi = mid
            break

        if value < 0:
            lo = mid
        else:
            hi = mid

    assert C.eval(lo) <= 0
    assert C.eval(hi) >= 0

    intervals[k] = (lo, hi)

    print(
        f"s_{k}: interval=({lo}, {hi}), "
        f"approx={float((lo+hi)/2):.12g}",
        flush=True
    )

print("STAGE 3: verify root ordering", flush=True)

all_ordered = True

for k in range(14, 22):
    lo1, hi1 = intervals[k]
    lo2, hi2 = intervals[k+1]

    if hi1 < lo2:
        print(f"PASS: s_{k} < s_{k+1}")
    elif hi2 < lo1:
        print(f"COUNTEREXAMPLE: s_{k} > s_{k+1}")
        all_ordered = False
    else:
        print(
            f"STOP-GATE: s_{k} and s_{k+1} "
            "not separated"
        )
        all_ordered = False

print("STAGE 4: conclusion", flush=True)

if all_ordered:
    print("CERTIFIED: all nine transition roots ordered")
    print("CERTIFIED: one x-sign change for every u>0")
    print("CERTIFIED: one positive discriminant root")
else:
    print("STOP-GATE: universal ordering not certified")
    print("UNIVERSAL DISCRIMINANT ROOT COUNT: OPEN")

print("DISCRIMINANT/B1 GLOBAL ROOT ORDERING: OPEN")
print("GENERAL-y POSITIVITY: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("R29AQ COMPLETE")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3),
    flush=True
)
