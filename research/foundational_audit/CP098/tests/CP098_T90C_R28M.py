
import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28M ===", flush=True)

source = Path(__file__).with_name("CP098_T90C_R28K2.py")

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

at = ns["at"]
pairs = ns["pairs"]

B = ns["ns"]["B"]
H = ns["ns"]["H"]
phi1 = ns["ns"]["phi1"]

assert len(pairs) == 11

print("STAGE 1: exact witness", flush=True)

assert ns["bracket_at"](B, H) == 0
print("PASS: D=0", flush=True)

print("STAGE 2: canonical derivative", flush=True)

# Hamiltonian vector field of phi1, evaluated at
# the witness only after differentiation.
directions = {}

for x, p in pairs:
    directions[x] = at(S.diff(phi1, p))
    directions[p] = -at(S.diff(phi1, x))

def directional_at(expr):
    total = S.Integer(0)

    for variable, coefficient in directions.items():
        if coefficient == 0:
            continue

        total += coefficient * at(S.diff(expr, variable))

    return S.simplify(total)

result = S.Integer(0)

for index, (x, p) in enumerate(pairs, 1):
    print(
        "Canonical pair",
        index,
        "of",
        len(pairs),
        flush=True
    )

    Bx = S.diff(B, x)
    Bp = S.diff(B, p)
    Hx = S.diff(H, x)
    Hp = S.diff(H, p)

    contribution = (
        directional_at(Bx) * at(Hp)
        + at(Bx) * directional_at(Hp)
        - directional_at(Bp) * at(Hx)
        - at(Bp) * directional_at(Hx)
    )

    result += S.cancel(contribution)

result = S.simplify(result)

print("K_D1 =", result, flush=True)
print("R28M FIRST COEFFICIENT COMPLETE", flush=True)
print("K_D2: OPEN", flush=True)
print("F_D: OPEN", flush=True)
print("Ddot preservation: OPEN", flush=True)
print("Full Dirac closure: OPEN", flush=True)
print(
    "ELAPSED SECONDS =",
    round(perf_counter() - start, 3),
    flush=True
)
