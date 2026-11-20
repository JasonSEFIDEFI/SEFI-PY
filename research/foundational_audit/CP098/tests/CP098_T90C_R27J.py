import ast
import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27J ===", flush=True)

source = Path(
    "research/foundational_audit/CP098/tests/"
    "CP098_T90C_R27I2.py"
)

tree = ast.parse(source.read_text())

# Reuse the model and witness definitions only.
# Stop before R27I2's computational stages.
cut = None
for i, node in enumerate(tree.body):
    if isinstance(node, ast.Expr):
        call = node.value
        if isinstance(call, ast.Call):
            if isinstance(call.func, ast.Name):
                if call.func.id == "print":
                    if call.args and isinstance(
                        call.args[0], ast.Constant
                    ):
                        if str(call.args[0].value).startswith(
                            "STAGE 1:"
                        ):
                            cut = i
                            break

if cut is None:
    raise RuntimeError("Cannot locate R27I2 stage boundary")

prefix = ast.Module(body=tree.body[:cut], type_ignores=[])
namespace = {"__name__": "__r27j_model__"}
exec(compile(prefix, str(source), "exec"), namespace)

u = namespace["u"]
q = namespace["q"]
P = namespace["P"]
p1 = namespace["p1"]
B = namespace["B"]
phi1 = namespace["phi1"]
phi2 = namespace["phi2"]
A = namespace["A"]
J = namespace["J"]
W = namespace["W"]
v = namespace["v"]
cases = namespace["cases"]
ev = namespace["ev"]

print("STAGE 1: canonical directional derivatives", flush=True)

# Canonical pairs (u,p1), (q,P).
# Since phi1=u.P and phi2=q.P:
#
# {B,phi1} = u.dB/dq - P.dB/dp1
# {B,phi2} = q.dB/dq - P.dB/dP

K1_direct = sum(
    u[i]*S.diff(B,q[i])
    -P[i]*S.diff(B,p1[i])
    for i in range(4)
)

K2_direct = sum(
    q[i]*S.diff(B,q[i])
    -P[i]*S.diff(B,P[i])
    for i in range(4)
)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: verify singular witnesses", flush=True)

for label, data in cases:
    for name, expr in [
        ("A", A),
        ("B", B),
        ("phi1", phi1),
        ("phi2", phi2),
    ]:
        value = ev(expr, data)
        if value != 0:
            raise AssertionError(
                f"{label}: {name}={value}"
            )

    print("PASS:", label, flush=True)

print("STAGE 2 PASS", flush=True)

print("STAGE 3: compare independent coefficients", flush=True)

for label, data in cases:
    k1 = ev(K1_direct, data)
    k2 = ev(K2_direct, data)
    predicted = ev(2*W*J/v**2, data)

    print(label, flush=True)
    print("  {B,phi1} =", k1, flush=True)
    print("  {B,phi2} =", k2, flush=True)
    print("  2WJ/v² =", predicted, flush=True)

    if k1 != 0:
        raise AssertionError(
            f"{label}: K1 should vanish, got {k1}"
        )

    if k2 != predicted:
        raise AssertionError(
            f"{label}: K2 mismatch: {k2} vs {predicted}"
        )

    print("  PASS: multiplier coefficients", flush=True)

print("STAGE 3 PASS", flush=True)

print()
print("R27J COMPLETE")
print("4D CANONICAL MULTIPLIER AUDIT: PASS")
print("SINGULAR WITNESS COEFFICIENTS: VERIFIED")
print("GENERAL SINGULAR CLASSIFICATION: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3)
)
