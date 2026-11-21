import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29AM ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29AL.py").read_text()

marker = '\nprint("STAGE 1: reconstruct B1"'
assert source.count(marker) == 1, "STOP: ambiguous source boundary"

ns = {
    "__name__": "r29am_definitions",
    "__file__": str(folder / "CP098_T90C_R29AL.py")
}
exec(
    compile(source[:source.index(marker)], "R29AL definitions", "exec"),
    ns
)

x, y, h, b = ns["x"], ns["y"], ns["h"], ns["b"]
ND, NT = ns["ND"], ns["NT"]
reduce_expr = ns["reduce_expr"]

q, u = S.symbols("q u")

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
    assert ex >= 4 and eh % 2 == 0
    terms[(ex - 4, ey, eh // 2)] = coeff

G = S.Poly.from_dict(terms, (x, y, q))
assert G.degree(q) == 3

Gq = S.Poly(G.as_expr(), q)
a0, a1, a2, a3 = [Gq.nth(k) for k in range(4)]

print("PASS: general cubic reconstructed", flush=True)

print("STAGE 2: construct exact discriminant", flush=True)
t0 = perf_counter()

Delta = S.Poly(
    a1**2*a2**2
    - 4*a3*a1**3
    - 4*a0*a2**3
    - 27*a3**2*a0**2
    + 18*a3*a2*a1*a0,
    x, y
)

assert not Delta.is_zero

print("DISCRIMINANT TERMS =", len(Delta.terms()))
print("DEGREES x,y =", Delta.degree(x), Delta.degree(y))
print("SECONDS =", round(perf_counter()-t0, 3))

print("STAGE 3: exact x-factor removal", flush=True)

min_x = min(ex for (ex, ey), c in Delta.terms())
print("MINIMUM x POWER =", min_x)

reduced = S.Poly.from_dict(
    {
        (ex-min_x, ey): c
        for (ex, ey), c in Delta.terms()
    },
    (x, y)
)

assert S.Poly(
    x**min_x * reduced.as_expr(), x, y
) == Delta

print("PASS: discriminant factorization")

print("STAGE 4: admissible-boundary shift", flush=True)
t0 = perf_counter()

shifted = S.Poly(
    reduced.as_expr().subs(y, u + S.Rational(31,32)),
    x, u
)

pos = sum(1 for c in shifted.coeffs() if c.is_positive is True)
neg = sum(1 for c in shifted.coeffs() if c.is_negative is True)

assert pos + neg == len(shifted.terms()), (
    "STOP: unresolved coefficient sign"
)

print("SHIFTED TERMS =", len(shifted.terms()))
print("SHIFTED DEGREES x,u =", shifted.degree(x), shifted.degree(u))
print("POSITIVE COEFFICIENTS =", pos)
print("NEGATIVE COEFFICIENTS =", neg)
print("SHIFT SECONDS =", round(perf_counter()-t0, 3))

print("STAGE 5: inspect x-coefficient sign groups", flush=True)

Px = S.Poly(shifted.as_expr(), x)

for k in range(Px.degree()+1):
    C = S.Poly(Px.nth(k), u)
    if C.is_zero:
        continue

    cp = sum(1 for c in C.coeffs() if c.is_positive is True)
    cn = sum(1 for c in C.coeffs() if c.is_negative is True)

    assert cp + cn == len(C.terms())

    if cp and not cn:
        status = "POSITIVE"
    elif cn and not cp:
        status = "NEGATIVE"
    else:
        status = "MIXED"

    print(
        f"x^{k}: degree_u={C.degree()} "
        f"positive={cp} negative={cn} "
        f"status={status}",
        flush=True
    )

print("R29AM COMPLETE")
print("GENERAL DISCRIMINANT STRUCTURE: VERIFIED")
print("DISCRIMINANT SIGN OVERLAP: OPEN")
print("GENERAL-y POSITIVITY: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =", round(perf_counter()-start, 3))
