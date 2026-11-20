import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28E ===", flush=True)

source = Path(__file__).with_name(
    "CP098_T90C_R27U.py"
)
assert source.is_file(), "R27U source missing"

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

w = S.symbols("w", real=True, nonzero=True)
rho, t = S.symbols("rho t", real=True)

u, q, P = ns["u"], ns["q"], ns["P"]
p0, p1 = ns["p0"], ns["p1"]

point = dict(ns["point"])

point[q[1]] = w
point[P[2]] = rho
point[p1[3]] = t

W0 = w**2
z0 = rho**2

IS0 = 12/W0 + S.Rational(448, 961)*z0

point[ns["IS"]] = IS0

IA0 = -S.Rational(3, 2)
C0 = S.Rational(31, 8)

M0 = 1 + (IA0**2 + IS0**2 + 1)/2

point[p0[0]] = W0*z0/(4*C0) - M0

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

print("STAGE 2: full canonical bracket",
      flush=True)

pairs = (
    list(zip(u, p1))
    + list(zip(q, P))
    + [
        (ns["IA"], ns["pA"]),
        (ns["IS"], ns["pS"]),
        (ns["IW"], ns["pW"])
    ]
)

assert len(pairs) == 11

parts = []

for i, (x, p) in enumerate(pairs, 1):
    Bx = at(S.diff(ns["B"], x))
    Bp = at(S.diff(ns["B"], p))
    Hx = at(S.diff(ns["H"], x))
    Hp = at(S.diff(ns["H"], p))

    term = S.cancel(Bx*Hp - Bp*Hx)
    parts.append(term)

    print(f"pair {i:02d} = {term}",
          flush=True)

D = S.cancel(sum(parts))

print("D(w,rho,t) =", S.factor(D),
      flush=True)

print("STAGE 2 PASS", flush=True)

print("STAGE 3: exact identity regression",
      flush=True)

baseline = (
    476656*W0**3*z0
    + 802816*W0**2*z0**2
    + 19393941*W0**2
    + 41330688*W0*z0
    + 531948096
)/(14776336*W0)

correction = S.Rational(8, 31)*W0*t**2

difference = S.cancel(
    D - baseline - correction
)

print("identity residual =", difference,
      flush=True)

assert difference == 0, (
    "Proposed t correction is not exact",
    difference
)

print("PASS: exact three-variable identity")

assert S.cancel(D.subs(t, 0)-baseline) == 0

print("PASS: R28C recovered at t=0")

checks = [
    (
        S.Integer(1),
        S.Integer(1),
        S.Integer(1),
        S.Rational(597765445, 14776336)
    ),
    (
        S.Integer(2),
        S.Integer(1),
        S.Integer(1),
        S.Rational(69496057, 3694084)
    ),
    (
        S.Rational(1, 2),
        S.Integer(2),
        S.Integer(1),
        S.Rational(9208989445, 59105344)
    ),
    (
        S.Integer(1),
        S.Integer(0),
        S.Integer(1),
        S.Rational(18635, 496)
    )
]

for wv, rv, tv, expected in checks:
    actual = S.cancel(D.subs({
        w: wv,
        rho: rv,
        t: tv
    }))
    assert actual == expected, (
        wv, rv, tv, actual, expected
    )

print("PASS: R28D witness regressions")
print("STAGE 3 PASS", flush=True)

print("STAGE 4: positivity certificate",
      flush=True)

# W=w²>0, z=rho²>=0, t²>=0.
# Every term in the R28C numerator is
# nonnegative, and its constant is positive.

W, z, tau = S.symbols(
    "W z tau",
    positive=True
)

N = (
    476656*W**3*z
    + 802816*W**2*z**2
    + 19393941*W**2
    + 41330688*W*z
    + 531948096
)

D_positive_form = (
    N/(14776336*W)
    + S.Rational(8, 31)*W*tau
)

assert S.cancel(
    D - D_positive_form.subs({
        W: w**2,
        z: rho**2,
        tau: t**2
    })
) == 0

assert all(
    coeff > 0
    for coeff in S.Poly(N, W, z).coeffs()
)

assert S.Integer(531948096) > 0

print("PASS: numerator strictly positive")
print("PASS: denominator strictly positive")
print("PASS: transverse correction nonnegative")
print("PASS: D>0 for all real rho,t and w!=0")
print("STAGE 4 PASS", flush=True)

print()
print("R28E COMPLETE")
print("THREE-VARIABLE IDENTITY: VERIFIED")
print("JOINT FAMILY POSITIVITY: PROVED")
print("GENERAL IA BRANCH: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3)
)
