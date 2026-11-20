import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29D ===", flush=True)

folder = Path(__file__).parent

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(folder / "CP098_T90C_R27U.py"))

u, q, P, p0, p1 = (
    ns[k] for k in ("u", "q", "P", "p0", "p1")
)

B = ns["B"]
IA, IS, IW, pA, pS, pW = (
    ns[k] for k in ("IA", "IS", "IW", "pA", "pS", "pW")
)

w, rho, z, b = S.symbols("w rho z b", real=True)

ISv = 12/w**2 + S.Rational(448, 961)*rho**2

Mv = 1 + (
    S.Rational(9, 4) + ISv**2 + 1 + b**2
)/2

point = {}

for vec, values in (
    (u, [2, 0, 0, 0]),
    (q, [0, w, 0, 0]),
    (P, [0, 0, rho, 0]),
    (p1, [0, 0, z, 0]),
    (p0, [
        2*w**2*rho**2/31 - Mv,
        -8*w*rho*z/31,
        0,
        0
    ])
):
    point.update({
        sym: S.sympify(val)
        for sym, val in zip(vec, values)
    })

point.update({
    IA: S.Rational(-3, 2),
    IS: ISv,
    IW: S.Integer(1),
    pA: S.Integer(0),
    pS: S.Integer(0),
    pW: b
})

print("STAGE 1: substitute general-w family", flush=True)

B_family = S.factor(S.cancel(B.subs(point)))

print("B FAMILY =", B_family, flush=True)

print("STAGE 2: verify B=0 identity", flush=True)

assert B_family == 0
assert S.simplify(S.diff(B_family, b)) == 0

print("PASS: B=0 for arbitrary b", flush=True)

print("STAGE 3: w=1 regression", flush=True)

expected = -S.Rational(32, 31)*rho*z

assert S.simplify(
    B_family.subs({w: 1, b: expected})
) == 0

print("PASS: previous restricted point satisfies B=0")
print("NOTE: B=0 does not determine b")
print("ORIGIN OF b RESTRICTION: OPEN")

print("R29D PASS")
print("GENERAL-w B=0 BRANCH: DERIVED")
print("GENERAL-w D=0: OPEN")
print("GENERAL-w Ddot: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3)
)
