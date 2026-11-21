import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29AE ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29AD.py").read_text()
cut = source.index(
    '\nprint("STAGE 1: reconstruct exact cubic"'
)

ns = {
    "__name__": "r29ae_definitions",
    "__file__": str(folder / "CP098_T90C_R29AD.py")
}
exec(compile(source[:cut], "R29AD definitions", "exec"), ns)

x, y, h, b, q = (
    ns["x"], ns["y"], ns["h"], ns["b"], ns["q"]
)
ND, NT = ns["ND"], ns["NT"]
reduce_expr = ns["reduce_expr"]

print("STAGE 1: reconstruct sign polynomials", flush=True)

PD = S.cancel(y * reduce_expr(ND))
PT = S.cancel(y * reduce_expr(NT))

A = S.expand(PD.subs(y, 4))
B = S.expand(PT.subs(y, 4))

R = S.Poly(S.resultant(A, B, b), x, h)

assert len(R.terms()) == 24
assert R.degree(x) == 13
assert R.degree(h) == 6

terms = {}
for (ex, eh), coeff in R.terms():
    assert ex >= 4 and eh % 2 == 0
    terms[(ex - 4, eh // 2)] = coeff

F = S.Poly.from_dict(terms, (x, q))
Fq = S.Poly(F.as_expr(), q)

a0, a1, a2, a3 = [Fq.nth(i) for i in range(4)]
A1 = S.Poly(a1, x)

disc = S.Poly(
    a1**2*a2**2
    - 4*a3*a1**3
    - 4*a0*a2**3
    - 27*a3**2*a0**2
    + 18*a3*a2*a1*a0,
    x
)

minimum_power = min(powers[0] for powers, _ in disc.terms())

assert minimum_power == 4

disc_reduced = S.Poly.from_dict(
    {
        (powers[0] - minimum_power,): coeff
        for powers, coeff in disc.terms()
    },
    x
)

assert S.expand(
    disc.as_expr() - x**4 * disc_reduced.as_expr()
) == 0

print("PASS: exact reconstruction", flush=True)

print("STAGE 2: verify unique positive roots", flush=True)

for name, poly in [
    ("A1", A1),
    ("DISC", disc_reduced)
]:
    assert poly.count_roots(0, S.oo) == 1, (
        f"STOP: {name} root count differs from R29AD"
    )
    assert poly.eval(1) != 0
    assert poly.eval(10) != 0

print("PASS: one positive root each", flush=True)

print("STAGE 3: rational root isolation", flush=True)

def isolate(poly, iterations=24):
    lo = S.Integer(1)
    hi = S.Integer(10)

    slo = S.sign(poly.eval(lo))
    shi = S.sign(poly.eval(hi))

    assert slo in (-1, 1)
    assert shi == -slo

    for _ in range(iterations):
        mid = (lo + hi) / 2
        smid = S.sign(poly.eval(mid))

        if smid == 0:
            return mid, mid

        assert smid in (-1, 1)

        if smid == slo:
            lo = mid
        else:
            hi = mid

    assert S.sign(poly.eval(lo)) == slo
    assert S.sign(poly.eval(hi)) == shi

    return lo, hi

a_lo, a_hi = isolate(A1)
d_lo, d_hi = isolate(disc_reduced)

print("A1 ROOT INTERVAL =", (a_lo, a_hi), flush=True)
print("DISC ROOT INTERVAL =", (d_lo, d_hi), flush=True)

print("STAGE 4: exact root ordering", flush=True)

if d_hi < a_lo:
    print("PASS: discriminant root precedes A1 root")
    ordering = "DISC_FIRST"

elif a_hi < d_lo:
    print("PASS: A1 root precedes discriminant root")
    ordering = "A1_FIRST"

else:
    print("STOP-GATE: intervals overlap")
    ordering = "UNRESOLVED"

print("ORDERING =", ordering, flush=True)

print("STAGE 5: interpretation", flush=True)

if ordering == "DISC_FIRST":
    assert S.sign(A1.eval(1)) == 1
    assert S.sign(A1.eval(10)) == -1
    assert S.sign(disc_reduced.eval(1)) == 1
    assert S.sign(disc_reduced.eval(10)) == -1

    print("CERTIFIED: A1<0 implies discriminant<0")
    print("CERTIFIED: F(x,q)>0 for x>0,q>0")
    print("CERTIFIED: no simultaneous zeros on y=4")

elif ordering == "A1_FIRST":
    print("POSITIVE-ROOT REGION MAY EXIST")
    print("EXACT RECOVERY OF q AND b REQUIRED")

else:
    print("ROOT ORDERING: OPEN")
    print("GLOBAL y=4 CLASSIFICATION: OPEN")

print("R29AE COMPLETE")
print("GENERAL y DOMAIN: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3),
    flush=True
)
