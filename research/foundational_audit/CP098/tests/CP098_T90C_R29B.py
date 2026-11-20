import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29B ===", flush=True)

folder = Path(__file__).parent

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(folder / "CP098_T90C_R27U.py"))

u, q, P, p0, p1 = (
    ns[k] for k in ("u", "q", "P", "p0", "p1")
)

B, H = ns["B"], ns["H"]

IA, IS, IW, pA, pS, pW = (
    ns[k] for k in ("IA", "IS", "IW", "pA", "pS", "pW")
)

pairs = (
    list(zip(u, p1))
    + list(zip(q, P))
    + [(IA, pA), (IS, pS), (IW, pW)]
)

rho, z = S.symbols("rho z", real=True)

b = -S.Rational(32, 31) * rho * z

ISv = 12 + S.Rational(448, 961) * rho**2

Mv = 1 + (
    S.Rational(9, 4)
    + ISv**2
    + 1
    + b**2
) / 2

point = {}

for vec, values in (
    (u, [2, 0, 0, 0]),
    (q, [0, 1, 0, 0]),
    (P, [0, 0, rho, 0]),
    (p1, [0, 0, z, 0]),
    (p0, [
        2 * rho**2 / 31 - Mv,
        -8 * rho * z / 31,
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

def at(expr):
    return S.cancel(expr.subs(point))

print("STAGE 1: Hamiltonian vector field", flush=True)

directions = {}

for x, p in pairs:
    directions[x] = at(S.diff(H, p))
    directions[p] = -at(S.diff(H, x))

active = [
    (var, val)
    for var, val in directions.items()
    if val != 0
]

print("Active directions:", len(active), flush=True)

def directional(expr):
    result = S.Integer(0)
    for var, val in active:
        result += val * at(S.diff(expr, var))
    return S.cancel(result)

print("STAGE 2: first canonical derivatives", flush=True)

first = [
    (
        S.diff(B, x),
        S.diff(B, p),
        S.diff(H, x),
        S.diff(H, p)
    )
    for x, p in pairs
]

print("PASS: first derivatives constructed", flush=True)

print("STAGE 3: second-bracket product rule", flush=True)

F = S.Integer(0)

for i, (Bx, Bp, Hx, Hp) in enumerate(first, 1):
    contribution = (
        directional(Bx) * at(Hp)
        + at(Bx) * directional(Hp)
        - directional(Bp) * at(Hx)
        - at(Bp) * directional(Hx)
    )

    F = S.cancel(F + contribution)

    print("F_D pair", i, "of 11", flush=True)

F = S.factor(F)

print("F_D UNREDUCED =", F, flush=True)

print("STAGE 4: impose D=0", flush=True)

Dbase = (
    802816 * rho**4
    + 41807344 * rho**2
    + 551342037
) / S.Integer(14776336)

z2 = S.factor(
    961 * Dbase / (8 * (32 * rho**2 - 31))
)

num, den = S.fraction(S.cancel(F))

relation = S.Poly(z**2 - z2, z)

num_red = S.rem(
    S.Poly(num, z), relation
).as_expr()

den_red = S.rem(
    S.Poly(den, z), relation
).as_expr()

num_red = S.factor(num_red)
den_red = S.factor(den_red)

print("F_D reduced numerator =", num_red, flush=True)
print("F_D reduced denominator =", den_red, flush=True)
print("D=0: z^2 =", z2, flush=True)

print("STAGE 5: second-witness regression", flush=True)

rho2 = S.Integer(2)

z2_second = S.factor(z2.subs(rho, rho2))
z_second = S.sqrt(z2_second)

expected = (
    -S.Rational(
        645017186759,
        33402088538116
    ) * S.sqrt(S.Integer(141894794986))
)

actual = S.simplify(
    F.subs({rho: rho2, z: z_second})
)

print("F_D second witness =", actual, flush=True)

if S.simplify(actual - expected) != 0:
    print("STOP: R28S regression mismatch", flush=True)
    raise AssertionError("F_D does not match R28S")

print("PASS: R28S exact regression", flush=True)
print("R29B PASS", flush=True)
print("COMBINED Ddot: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")

print(
    "ELAPSED SECONDS =",
    round(perf_counter() - start, 3),
    flush=True
)
