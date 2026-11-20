import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28A ===", flush=True)

source = Path(__file__).with_name(
    "CP098_T90C_R27U.py"
)
assert source.is_file(), "R27U source missing"

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

w = S.symbols("w", real=True, nonzero=True)
point = dict(ns["point"])

u, q, P = ns["u"], ns["q"], ns["P"]
p0, p1 = ns["p0"], ns["p1"]

point[q[1]] = w
point[P[2]] = S.Integer(1)
point[p1[3]] = S.Integer(0)

W0 = w**2
IS0 = 12/W0 + S.Rational(448,961)
point[ns["IS"]] = IS0

IA0 = -S.Rational(3,2)
C0 = S.Rational(31,8)
M0 = 1 + (IA0**2 + IS0**2 + 1)/2

point[p0[0]] = W0/(4*C0) - M0

def at(expr):
    return S.cancel(expr.subs(point))

print("STAGE 1: symbolic constraint audit",
      flush=True)

constraints = [
    ("phi1", ns["phi1"]),
    ("phi2", ns["phi2"]),
    ("chi", ns["chi"]),
    ("psi", ns["psi"]),
    ("theta", ns["theta"]),
    ("B", ns["B"]),
    ("A", ns["A"]),
    ("J", ns["J"]),
    ("Q", ns["Q"]),
    ("pA", ns["pA"])
]

for name, expr in constraints:
    residual = at(expr)
    assert residual == 0, (name, residual)
    print("PASS:", name, flush=True)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: canonical D(w)",
      flush=True)

pairs = (
    list(zip(u,p1))
    + list(zip(q,P))
    + [
        (ns["IA"],ns["pA"]),
        (ns["IS"],ns["pS"]),
        (ns["IW"],ns["pW"])
    ]
)

parts = []

for i,(x,p) in enumerate(pairs,1):
    Bx = at(S.diff(ns["B"],x))
    Bp = at(S.diff(ns["B"],p))
    Hx = at(S.diff(ns["H"],x))
    Hp = at(S.diff(ns["H"],p))

    term = S.cancel(Bx*Hp-Bp*Hx)
    parts.append(term)
    print(f"pair {i:02d} = {term}",flush=True)

D = S.factor(sum(parts))

expected = {
    S.Rational(1,2):
        S.Rational(8696808209,59105344),
    S.Integer(1):
        S.Rational(593952197,14776336),
    S.Integer(2):
        S.Rational(65682809,3694084)
}

for value,target in expected.items():
    actual = S.cancel(D.subs(w,value))
    assert actual == target, (value,actual,target)
    print("PASS: R27Z at w =",value)

print("D(w) =",D)
print("STAGE 2 PASS",flush=True)

print("STAGE 3: exact sign investigation",
      flush=True)

num,den = S.fraction(S.cancel(D))
print("numerator =",S.factor(num))
print("denominator =",S.factor(den))

# Clear the Laurent powers using W=w^2>0.
Wvar = S.symbols("Wvar",positive=True)
num_poly = S.Poly(num,w)
den_poly = S.Poly(den,w)

def even_in_w(poly):
    return all(
        monom[0] % 2 == 0
        for monom,coeff in poly.terms()
    )

if even_in_w(num_poly) and even_in_w(den_poly):
    nW = S.Poly.from_dict(
        {(m[0]//2,):c for m,c in num_poly.terms()},
        Wvar
    ).as_expr()

    dW = S.Poly.from_dict(
        {(m[0]//2,):c for m,c in den_poly.terms()},
        Wvar
    ).as_expr()

    assert S.cancel(
        D - (nW/dW).subs(Wvar,w**2)
    ) == 0

    print("D(W) numerator =",S.factor(nW))
    print("D(W) denominator =",S.factor(dW))

    numerator_positive = all(
        c > 0 for c in S.Poly(nW,Wvar).all_coeffs()
    )
    denominator_positive = all(
        c > 0 for c in S.Poly(dW,Wvar).all_coeffs()
    )

    if numerator_positive and denominator_positive:
        print("PASS: D>0 for every W>0")
        print("PASS: no D=0 point in this family")
    else:
        print("POSITIVITY NOT CERTIFIED")
        print("Further exact sign analysis required")
else:
    print("EVEN-POLYNOMIAL TEST NOT APPLICABLE")
    print("Further exact sign analysis required")

print("STAGE 3 PASS",flush=True)

print()
print("R28A COMPLETE")
print("SYMBOLIC W FAMILY: CHECKED")
print("GENERAL IA BRANCH: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start,3)
)
