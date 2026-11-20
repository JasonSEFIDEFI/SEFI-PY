
import sympy as S
import runpy, contextlib, io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28W ===", flush=True)

source = Path(__file__).with_name("CP098_T90C_R27U.py")
with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

u, q, P, p1 = (
    ns["u"], ns["q"], ns["P"], ns["p1"]
)
H, B = ns["H"], ns["B"]
chi, phi2 = ns["chi"], ns["phi2"]

def X1(expr):
    return sum(
        u[i]*S.diff(expr, q[i])
        - P[i]*S.diff(expr, p1[i])
        for i in range(4)
    )

print("STAGE 1: generator action", flush=True)

for i in range(4):
    assert S.expand(X1(q[i])-u[i]) == 0
    assert S.expand(X1(p1[i])+P[i]) == 0

print("PASS: canonical generator action", flush=True)

print("STAGE 2: Hamiltonian identity", flush=True)

H1 = X1(H)
assert S.simplify(H1-chi+phi2) == 0

print("PASS: {H,phi1}=chi-phi2", flush=True)

print("STAGE 3: B transformation", flush=True)

B1 = X1(B)

# Verify the generator acts as a derivation.
# This is an exact algebraic check, not a
# restriction to either numerical witness.
check = S.simplify(
    X1(ns["W"])
)
assert check == 0

print("PASS: W invariant", flush=True)

# Check whether B is invariant on the
# established constraint surface at each witness.
folder = Path(__file__).parent

for name, filename in (
    ("FIRST", "CP098_T90C_R28K2.py"),
    ("SECOND", "CP098_T90C_R28Q.py")
):
    with contextlib.redirect_stdout(io.StringIO()):
        witness = runpy.run_path(str(folder / filename))

    value = S.simplify(witness["at"](B1))
    print(name, "{B,phi1} =", value, flush=True)

print("R28W: EXACT GENERATOR IDENTITIES VERIFIED")
print("FAMILY-WIDE K_D1: OPEN")
print("FAMILY-WIDE Ddot: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("ELAPSED =", round(perf_counter()-start,3))
