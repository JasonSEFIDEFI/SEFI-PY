import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29AO ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29AN.py").read_text()

marker = '\nprint("STAGE 1: reconstruct exact polynomials"'
assert source.count(marker) == 1, (
    "STOP: ambiguous source boundary"
)

ns = {
    "__name__": "r29ao_definitions",
    "__file__": str(folder / "CP098_T90C_R29AN.py")
}

exec(
    compile(source[:source.index(marker)],
            "R29AN definitions", "exec"),
    ns
)

x, y, h, b = ns["x"], ns["y"], ns["h"], ns["b"]
ND, NT = ns["ND"], ns["NT"]
reduce_expr = ns["reduce_expr"]

q, u = S.symbols("q u")

print("STAGE 1: reconstruct cubic coefficients", flush=True)

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

print("PASS: exact cubic reconstructed", flush=True)

print("STAGE 2: verify positive coefficient families", flush=True)

for k, ak in ((0, a0), (2, a2), (3, a3)):
    shifted = S.Poly(
        ak.subs(y, u + S.Rational(31, 32)),
        x, u
    )

    assert not shifted.is_zero
    assert all(
        c.is_positive is True
        for c in shifted.coeffs()
    ), f"STOP: B{k} positivity not certified"

    print(f"PASS: B{k}>0 for x>0,u>0", flush=True)

print("STAGE 3: exact discriminant identity", flush=True)

Delta = (
    a1**2*a2**2
    - 4*a3*a1**3
    - 4*a0*a2**3
    - 27*a3**2*a0**2
    + 18*a3*a2*a1*a0
)

Delta_at_B1_zero = -4*a0*a2**3 - 27*a3**2*a0**2

# Verify the identity using independent symbols first.
A0, A1, A2, A3 = S.symbols("A0 A1 A2 A3")

formal_delta = (
    A1**2*A2**2
    - 4*A3*A1**3
    - 4*A0*A2**3
    - 27*A3**2*A0**2
    + 18*A3*A2*A1*A0
)

formal_at_zero = -4*A0*A2**3 - 27*A3**2*A0**2

formal_quotient = (
    A1*A2**2
    - 4*A3*A1**2
    + 18*A3*A2*A0
)

assert S.Poly(
    formal_delta - formal_at_zero
    - A1*formal_quotient,
    A0, A1, A2, A3
).is_zero

# Verify the same identity with the actual
# reconstructed polynomial coefficients.
quotient = (
    a1*a2**2
    - 4*a3*a1**2
    + 18*a3*a2*a0
)

assert S.Poly(
    Delta - Delta_at_B1_zero - a1*quotient,
    x, y
).is_zero

print("PASS: exact discriminant identity", flush=True)

print("STAGE 4: boundary conclusion", flush=True)

print("CERTIFIED: Delta<0 whenever B1=0")
print("CERTIFIED: B1=0 and Delta=0 cannot intersect")
print("ADDITIONAL DISCRIMINANT ROOTS: NOT EXCLUDED")
print("GENERAL-y POSITIVITY: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")

print("R29AO COMPLETE")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3),
    flush=True
)
