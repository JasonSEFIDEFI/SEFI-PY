import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29AP ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29AO.py").read_text()

marker = '\nprint("STAGE 1: reconstruct cubic coefficients"'
assert source.count(marker) == 1, "STOP: source boundary ambiguous"

ns = {
    "__name__": "r29ap_definitions",
    "__file__": str(folder / "CP098_T90C_R29AO.py")
}

exec(
    compile(source[:source.index(marker)],
            "R29AO definitions", "exec"),
    ns
)

x, y, h, b = ns["x"], ns["y"], ns["h"], ns["b"]
ND, NT = ns["ND"], ns["NT"]
reduce_expr = ns["reduce_expr"]

q, u = S.symbols("q u")

print("STAGE 1: reconstruct exact discriminant", flush=True)

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
    {(ex - min_x, ey): c
     for (ex, ey), c in Delta.terms()},
    (x, y)
)

shifted = S.Poly(
    D.as_expr().subs(y, u + S.Rational(31, 32)),
    x, u
)

assert len(shifted.terms()) == 833

print("PASS: exact shifted discriminant", flush=True)

print("STAGE 2: analyze mixed x-coefficients", flush=True)

Px = S.Poly(shifted.as_expr(), x)
unresolved = []

for k in range(14, 23):
    C = S.Poly(Px.nth(k), u)
    assert C.degree() >= 0
    assert C.eval(0) != 0
    assert C.LC() != 0

    signs = []
    for c in C.all_coeffs():
        if c > 0:
            signs.append(1)
        elif c < 0:
            signs.append(-1)
        else:
            signs.append(0)

    nonzero = [s for s in signs if s != 0]
    changes = sum(
        nonzero[i] != nonzero[i - 1]
        for i in range(1, len(nonzero))
    )

    at_zero = S.sign(C.eval(0))
    at_infinity = S.sign(C.LC())

    print(
        f"x^{k}: degree_u={C.degree()}, "
        f"Descartes_changes={changes}, "
        f"sign_at_0={at_zero}, "
        f"sign_at_infinity={at_infinity}",
        flush=True
    )

    if changes == 0:
        print("  CERTIFIED: no positive u roots")
    elif changes == 1 and at_zero != at_infinity:
        print("  CERTIFIED: exactly one positive u root")
    else:
        print("  OPEN: exact positive u root count")
        unresolved.append(k)

print("STAGE 3: decision gate", flush=True)

if unresolved:
    print("UNRESOLVED COEFFICIENTS =", unresolved)
    print("STOP-GATE: further exact root analysis required")
else:
    print("PASS: every mixed coefficient root count certified")

print("R29AP COMPLETE")
print("DISCRIMINANT SIGN OVERLAP: OPEN")
print("GENERAL-y POSITIVITY: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter() - start, 3),
    flush=True
)
