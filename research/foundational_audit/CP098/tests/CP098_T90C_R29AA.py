import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29AA ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29Z.py").read_text()
cut = source.index('print("STAGE 1: reconstruct reduced resultant"')

ns = {
    "__name__": "r29aa_definitions",
    "__file__": str(folder / "CP098_T90C_R29Z.py")
}

exec(compile(source[:cut], "R29Z definitions", "exec"), ns)

x, y, h, b, q = (
    ns["x"], ns["y"], ns["h"], ns["b"], ns["q"]
)

ND, NT = ns["ND"], ns["NT"]
reduce_expr = ns["reduce_expr"]

print("STAGE 1: reconstruct exact certificate", flush=True)

PD = S.cancel(y * reduce_expr(ND))
PT = S.cancel(y * reduce_expr(NT))

A = S.expand(PD.subs(y, 4))
B = S.expand(PT.subs(y, 4))

R = S.Poly(S.resultant(A, B, b), x, h)

assert len(R.terms()) == 24
assert R.degree(x) == 13
assert R.degree(h) == 6

terms = {}
for (ex, eh), coefficient in R.terms():
    assert ex >= 4 and eh % 2 == 0
    terms[(ex - 4, eh // 2)] = coefficient

F = S.Poly.from_dict(terms, (x, q))
Fq = S.Poly(F.as_expr(), q)

A0 = S.Poly(Fq.nth(0), x)
A1 = S.Poly(Fq.nth(1), x)
A2 = S.Poly(Fq.nth(2), x)
A3 = S.Poly(Fq.nth(3), x)

C = S.Poly(
    4*A0.as_expr()*A2.as_expr() - A1.as_expr()**2,
    x
)

assert C.degree() == 16
assert len(C.terms()) == 17

print("PASS: R29Z certificate reconstructed", flush=True)

print("STAGE 2: exact rational sign evaluations", flush=True)

test_points = [
    S.Rational(1, 100),
    S.Rational(1, 10),
    S.Integer(1),
    S.Integer(10),
    S.Integer(100)
]

negative_witnesses = []

for value in test_points:
    result = C.eval(value)
    assert result.is_Rational

    sign = S.sign(result)
    assert sign in (S.Integer(-1), S.Integer(0), S.Integer(1))

    print(
        "x =", value,
        "C(x) sign =", sign,
        flush=True
    )

    if sign == -1:
        negative_witnesses.append(value)

print("STAGE 3: interpret certificate", flush=True)

if negative_witnesses:
    witness = negative_witnesses[0]
    assert C.eval(witness) < 0

    print(
        "PROVED: C(x)>0 is false globally",
        flush=True
    )
    print(
        "EXACT NEGATIVE WITNESS: x =", witness,
        flush=True
    )
    print(
        "FULL CUBIC POSITIVITY: STILL OPEN",
        flush=True
    )
else:
    print(
        "NO NEGATIVE WITNESS IN TEST SET",
        flush=True
    )
    print(
        "GLOBAL CERTIFICATE SIGN: OPEN",
        flush=True
    )

print("R29AA COMPLETE")
print("y=4 RESULTANT ZERO CLASSIFICATION: OPEN")
print("GENERAL y DOMAIN: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter() - start, 3),
    flush=True
)
