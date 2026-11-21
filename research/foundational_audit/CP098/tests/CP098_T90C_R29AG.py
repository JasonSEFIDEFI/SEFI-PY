import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29AG ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29AF.py").read_text()

marker = '\nprint("STAGE 1: reconstruct reduced equations"'
assert source.count(marker) == 1, (
    "STOP: R29AF source boundary ambiguous"
)

cut = source.index(marker)

ns = {
    "__name__": "r29ag_definitions",
    "__file__": str(folder / "CP098_T90C_R29AF.py")
}

exec(
    compile(source[:cut], "R29AF definitions", "exec"),
    ns
)

x = ns["x"]
y = ns["y"]
h = ns["h"]
b = ns["b"]
ND = ns["ND"]
NT = ns["NT"]
reduce_expr = ns["reduce_expr"]

q = S.Symbol("q", positive=True)

print("STAGE 1: construct exact equations", flush=True)

PD = S.Poly(S.cancel(y * reduce_expr(ND)), b)
PT = S.Poly(S.cancel(y * reduce_expr(NT)), b)

assert PD.degree() == 2
assert PT.degree() == 3

print("PASS: quadratic/cubic structure", flush=True)

print("STAGE 2: symbolic resultant", flush=True)

t0 = perf_counter()

# Subresultant-based elimination with symbolic
# parameters x,y,h; no numerical substitution.
R_expr = S.resultant(
    PD.as_expr(),
    PT.as_expr(),
    b
)

R = S.Poly(R_expr, x, y, h)

assert not R.is_zero

print(
    "PASS: general-y resultant constructed",
    flush=True
)
print(
    "RESULTANT TERMS =", len(R.terms()),
    flush=True
)
print(
    "DEGREES x,y,h =",
    R.degree(x), R.degree(y), R.degree(h),
    flush=True
)
print(
    "RESULTANT SECONDS =",
    round(perf_counter()-t0, 3),
    flush=True
)

print("STAGE 3: exact structural inspection", flush=True)

min_x = min(powers[0] for powers, _ in R.terms())

odd_h = [
    (powers, coefficient)
    for powers, coefficient in R.terms()
    if powers[2] % 2 != 0
]

print("MINIMUM x POWER =", min_x, flush=True)
print("ODD h TERMS =", len(odd_h), flush=True)

if min_x < 4 or odd_h:
    print("STOP: proposed universal reduction unavailable")
    print("GENERAL-y ZERO CLASSIFICATION: OPEN")
    raise SystemExit(0)

print("STAGE 4: exact reduction", flush=True)

reduced_terms = {}

for (ex, ey, eh), coefficient in R.terms():
    assert ex >= 4
    assert eh % 2 == 0

    key = (ex - 4, ey, eh // 2)
    reduced_terms[key] = coefficient

G = S.Poly.from_dict(reduced_terms, (x, y, q))

assert S.expand(
    R.as_expr()
    - x**4 * G.as_expr().subs(q, h**2)
) == 0

print("PASS: R = x^4 G(x,y,h^2)", flush=True)

print("G TERMS =", len(G.terms()), flush=True)
print(
    "G DEGREES x,y,q =",
    G.degree(x), G.degree(y), G.degree(q),
    flush=True
)

print("STAGE 5: regression against R29AF", flush=True)

for value in (1, 2, 4, 8):
    # Recompute each slice independently.
    A = S.Poly(PD.as_expr().subs(y, value), b)
    B = S.Poly(PT.as_expr().subs(y, value), b)

    slice_resultant = S.Poly(
        S.resultant(A.as_expr(), B.as_expr(), b),
        x, h
    )

    reconstructed = S.Poly(
        x**4 * G.as_expr().subs({
            y: value,
            q: h**2
        }),
        x, h
    )

    assert slice_resultant == reconstructed, (
        f"STOP: regression mismatch at y={value}"
    )

    print(
        f"PASS: exact y={value} regression",
        flush=True
    )

print("R29AG PASS")
print("GENERAL-y STRUCTURAL REDUCTION: VERIFIED")
print("GENERAL-y POSITIVITY: OPEN")
print("GENERAL-y SIMULTANEOUS ZERO CLASSIFICATION: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3),
    flush=True
)
