
import sympy as S
import runpy, contextlib, io
from pathlib import Path

print("=== CP098-T90C-R28R ===", flush=True)

source = Path(__file__).with_name("CP098_T90C_R28Q.py")

with contextlib.redirect_stdout(io.StringIO()):
    witness = runpy.run_path(str(source))

base = witness["base"]
at = witness["at"]
bracket_at = witness["bracket_at"]

H = base["H"]
phi1 = base["phi1"]
phi2 = base["phi2"]
Q = base["Q"]
B = base["B"]

print("STAGE 1: verify second witness", flush=True)

assert at(Q) == 0
assert bracket_at(B, H) == 0
print("PASS: Q=0 and D=0", flush=True)

print("STAGE 2: Q tangency", flush=True)

F = bracket_at(Q, H)
K1 = bracket_at(Q, phi1)
K2 = bracket_at(Q, phi2)

print("F_Q =", F, flush=True)
print("K_Q1 =", K1, flush=True)
print("K_Q2 =", K2, flush=True)

assert K1 == 0, "STOP: lambda1 enters Q tangency"
assert K2 != 0, "STOP: lambda2 not determined"

lam2 = S.simplify(-F/K2)

assert S.simplify(F + lam2*K2) == 0

print("lambda2 =", lam2, flush=True)
print("PASS: Qdot=0 fixes lambda2", flush=True)

print("STAGE 3: earlier tangencies", flush=True)

constraints = (
    "phi1", "phi2", "chi", "psi",
    "theta", "B", "A", "J", "pA"
)

for name in constraints:
    g = base[name]

    f = bracket_at(g, H)
    k1 = bracket_at(g, phi1)
    k2 = bracket_at(g, phi2)

    total = S.simplify(f + lam2*k2)

    print(name, "tangency =", total, flush=True)

    assert k1 == 0, (
        "STOP: lambda1 enters tangency of " + name
    )
    assert total == 0, (
        "STOP: tangency fails for " + name
    )

print("R28R PASS: ALL EARLIER TANGENCIES")
print("SECOND WITNESS lambda2 =", lam2)
print("Ddot PRESERVATION: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
