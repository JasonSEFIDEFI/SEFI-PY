import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29AF ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29AE.py").read_text()

# Anchor to the actual statement, not its quoted
# occurrence inside another source.index call.
marker = '\nprint("STAGE 1: reconstruct sign polynomials"'
assert source.count(marker) == 1
cut = source.index(marker)

ns = {
    "__name__": "r29af_definitions",
    "__file__": str(folder / "CP098_T90C_R29AE.py")
}
exec(compile(source[:cut], "R29AE definitions", "exec"), ns)

x, y, h, b = ns["x"], ns["y"], ns["h"], ns["b"]
ND, NT = ns["ND"], ns["NT"]
reduce_expr = ns["reduce_expr"]

print("STAGE 1: reconstruct reduced equations", flush=True)

PD = S.Poly(S.cancel(y * reduce_expr(ND)), b)
PT = S.Poly(S.cancel(y * reduce_expr(NT)), b)

assert PD.degree() == 2
assert PT.degree() == 3

print("PASS: reduced polynomial degrees", flush=True)

print("STAGE 2: inspect exact admissible slices", flush=True)

for yvalue in [1, 2, 4, 8]:
    t0 = perf_counter()
    print("SLICE y =", yvalue, flush=True)

    A = S.Poly(PD.as_expr().subs(y, yvalue), b)
    B = S.Poly(PT.as_expr().subs(y, yvalue), b)

    assert A.degree() == 2
    assert B.degree() == 3

    R = S.Poly(S.resultant(A, B, b), x, h)
    assert not R.is_zero

    terms = R.terms()
    min_x = min(powers[0] for powers, _ in terms)
    odd_h = sum(
        1 for powers, _ in terms
        if powers[1] % 2 != 0
    )

    print("  terms =", len(terms), flush=True)
    print("  degree x =", R.degree(x), flush=True)
    print("  degree h =", R.degree(h), flush=True)
    print("  minimum x power =", min_x, flush=True)
    print("  odd h terms =", odd_h, flush=True)

    assert odd_h == 0, (
        f"STOP: even-h structure failed at y={yvalue}"
    )

    assert min_x >= 0

    print(
        "  elapsed =",
        round(perf_counter() - t0, 3),
        flush=True
    )

print("R29AF COMPLETE")
print("EXACT STRUCTURE: VERIFIED ON FOUR SLICES")
print("OTHER SLICES: NOT CLASSIFIED")
print("GENERAL y DOMAIN: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter() - start, 3),
    flush=True
)
