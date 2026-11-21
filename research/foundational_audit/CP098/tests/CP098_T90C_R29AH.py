import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29AH ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29AG.py").read_text()

marker = '\nprint("STAGE 1: construct exact equations"'
assert source.count(marker) == 1, (
    "STOP: ambiguous source boundary"
)

ns = {
    "__name__": "r29ah_definitions",
    "__file__": str(folder / "CP098_T90C_R29AG.py")
}

exec(
    compile(
        source[:source.index(marker)],
        "R29AG definitions",
        "exec"
    ),
    ns
)

x, y, h, b = (
    ns["x"], ns["y"], ns["h"], ns["b"]
)

ND, NT = ns["ND"], ns["NT"]
reduce_expr = ns["reduce_expr"]

q = S.Symbol("q", positive=True)
u = S.Symbol("u", positive=True)

print("STAGE 1: reconstruct general resultant", flush=True)

PD = S.Poly(S.cancel(y * reduce_expr(ND)), b)
PT = S.Poly(S.cancel(y * reduce_expr(NT)), b)

assert PD.degree() == 2
assert PT.degree() == 3

R = S.Poly(
    S.resultant(PD.as_expr(), PT.as_expr(), b),
    x, y, h
)

assert len(R.terms()) == 63
assert R.degree(x) == 13
assert R.degree(y) == 11
assert R.degree(h) == 6

terms = {}

for (ex, ey, eh), coeff in R.terms():
    assert ex >= 4
    assert eh % 2 == 0
    terms[(ex - 4, ey, eh // 2)] = coeff

G = S.Poly.from_dict(terms, (x, y, q))

assert len(G.terms()) == 63
assert G.degree(x) == 9
assert G.degree(y) == 11
assert G.degree(q) == 3

print("PASS: R29AG reduction reconstructed", flush=True)

print("STAGE 2: shift admissible y boundary", flush=True)

H = S.Poly(
    G.as_expr().subs(y, u + S.Rational(31, 32)),
    x, u, q
)

assert S.expand(
    H.as_expr().subs(u, y - S.Rational(31, 32))
    - G.as_expr()
) == 0

print("PASS: exact boundary shift", flush=True)
print("SHIFTED TERMS =", len(H.terms()), flush=True)
print(
    "DEGREES x,u,q =",
    H.degree(x), H.degree(u), H.degree(q),
    flush=True
)

print("STAGE 3: coefficient sign inventory", flush=True)

positive = 0
negative = 0
unresolved = 0

for powers, coeff in H.terms():
    if coeff.is_positive is True:
        positive += 1
    elif coeff.is_negative is True:
        negative += 1
    else:
        unresolved += 1

print("POSITIVE TERMS =", positive, flush=True)
print("NEGATIVE TERMS =", negative, flush=True)
print("UNRESOLVED SIGNS =", unresolved, flush=True)

assert unresolved == 0, (
    "STOP: nonrational or unresolved coefficients"
)

print("STAGE 4: inspect q coefficients", flush=True)

for k in range(4):
    Bk = S.Poly(
        S.Poly(H.as_expr(), q).nth(k),
        x, u
    )

    pos = sum(
        1 for c in Bk.coeffs()
        if c.is_positive is True
    )
    neg = sum(
        1 for c in Bk.coeffs()
        if c.is_negative is True
    )

    assert pos + neg == len(Bk.terms())

    print(
        f"B{k}: terms={len(Bk.terms())}, "
        f"positive={pos}, negative={neg}",
        flush=True
    )

    if neg:
        print(
            f"B{k} NEGATIVE MONOMIAL EXPONENTS:",
            [
                powers
                for powers, coeff in Bk.terms()
                if coeff.is_negative is True
            ],
            flush=True
        )

print("STAGE 5: conclusion", flush=True)

if negative == 0:
    print(
        "CERTIFIED: G(x,y,q)>0 throughout "
        "x>0, y>31/32, q>0",
        flush=True
    )
    print(
        "CERTIFIED: no simultaneous zeros "
        "in the reduced admissible family",
        flush=True
    )
else:
    print(
        "STOP-GATE: direct shifted-coefficient "
        "positivity unavailable",
        flush=True
    )
    print(
        "GENERAL-y POSITIVITY: OPEN",
        flush=True
    )

print("R29AH COMPLETE")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter() - start, 3),
    flush=True
)
