
import sympy as S
import runpy, contextlib, io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28S ===", flush=True)

source = Path(__file__).with_name("CP098_T90C_R28Q.py")

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

base = ns["base"]
at = ns["at"]
bracket_at = ns["bracket_at"]

B = base["B"]
H = base["H"]
phi1 = base["phi1"]
phi2 = base["phi2"]

u, q, P = base["u"], base["q"], base["P"]
p0, p1 = base["p0"], base["p1"]

pairs = (
    list(zip(u, p1))
    + list(zip(q, P))
    + [
        (base["IA"], base["pA"]),
        (base["IS"], base["pS"]),
        (base["IW"], base["pW"])
    ]
)

assert len(pairs) == 11

print("STAGE 1: witness regression", flush=True)
assert bracket_at(B, H) == 0
print("PASS: D=0", flush=True)

# Compute exact first derivatives once.
print("STAGE 2: first derivatives", flush=True)

first = []
for x, p in pairs:
    first.append((
        S.diff(B, x),
        S.diff(B, p),
        S.diff(H, x),
        S.diff(H, p)
    ))

def coefficient(generator, label):
    print("STAGE:", label, flush=True)

    directions = {}
    for x, p in pairs:
        directions[x] = at(S.diff(generator, p))
        directions[p] = -at(S.diff(generator, x))

    active = [
        (variable, value)
        for variable, value in directions.items()
        if value != 0
    ]

    print("Active directions:", len(active), flush=True)

    def directional_at(expr):
        total = S.Integer(0)
        for variable, value in active:
            total += value * at(S.diff(expr, variable))
        return S.simplify(total)

    result = S.Integer(0)

    for i, (Bx, Bp, Hx, Hp) in enumerate(first, 1):
        print(label, "pair", i, "of 11", flush=True)

        term = (
            directional_at(Bx)*at(Hp)
            + at(Bx)*directional_at(Hp)
            - directional_at(Bp)*at(Hx)
            - at(Bp)*directional_at(Hx)
        )

        result += S.cancel(term)

    result = S.simplify(result)
    print(label, "=", result, flush=True)
    return result

K1 = coefficient(phi1, "K_D1")
K2 = coefficient(phi2, "K_D2")
F = coefficient(H, "F_D")

lam2 = (
    S.Rational(14,3201091)
    * S.sqrt(S.Integer(141894794986))
)

residual_without_lam1 = S.simplify(F + lam2*K2)

print("STAGE 3: preservation equation", flush=True)
print("K_D1 =", K1, flush=True)
print("K_D2 =", K2, flush=True)
print("F_D =", F, flush=True)
print(
    "Ddot without lambda1 =",
    residual_without_lam1,
    flush=True
)

if K1 == 0:
    if residual_without_lam1 == 0:
        print("PASS: Ddot=0 for arbitrary lambda1")
    else:
        print("OBSTRUCTION: Ddot != 0")
        print("lambda1 cannot cancel residual")
else:
    lam1 = S.simplify(-residual_without_lam1/K1)
    print("Ddot fixes lambda1 =", lam1)
    assert S.simplify(
        F + lam1*K1 + lam2*K2
    ) == 0
    print("PASS: Ddot=0 with determined lambda1")

print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3),
    flush=True
)
