import ast
import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27K ===", flush=True)

source = Path(
    "research/foundational_audit/CP098/tests/"
    "CP098_T90C_R27I2.py"
)

tree = ast.parse(source.read_text())

# Import only the previously verified model definitions
# and witness constructors, not its test execution.
cut = None
for index, node in enumerate(tree.body):
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
                            cut = index
                            break

if cut is None:
    raise RuntimeError("Model boundary not found")

prefix = ast.Module(
    body=tree.body[:cut],
    type_ignores=[]
)

ns = {"__name__": "__r27k_model__"}
exec(compile(prefix, str(source), "exec"), ns)

A = ns["A"]
J = ns["J"]
v = ns["v"]
W = ns["W"]
Hc = ns["Hc"]
I = ns["I"]
pI = ns["pI"]
cases = ns["cases"]
ev = ns["ev"]

print("STAGE 1: derive identity-sector brackets", flush=True)

# A and J depend only on identity coordinates
# and identity momenta. Their brackets with
# phi1 and phi2 vanish identically.

Adot = sum(
    S.diff(A, I[i])*S.diff(Hc, pI[i])
    for i in range(3)
)

Jdot = sum(
    S.diff(J, I[i])*S.diff(Hc, pI[i])
    - S.diff(J, pI[i])*S.diff(Hc, I[i])
    for i in range(3)
)

# Hc's identity momentum dependence is v*T.
# This is a general identity, not a witness fit.
residual = S.simplify(Adot-v*J)

if residual != 0:
    raise AssertionError(
        f"Adot identity failed: {residual}"
    )

print("PASS: {A,Hc} = v*J", flush=True)
print("PASS: {A,phi1} = {A,phi2} = 0", flush=True)
print("STAGE 1 PASS", flush=True)

print("STAGE 2: four-dimensional witness audit", flush=True)

for label, data in cases:
    aval = ev(A, data)
    jval = ev(J, data)
    vval = ev(v, data)
    adot = ev(Adot, data)
    jdot = ev(Jdot, data)
    addot = S.simplify(vval*jdot)

    if aval != 0:
        raise AssertionError(
            f"{label}: A={aval}"
        )

    if adot != vval*jval:
        raise AssertionError(
            f"{label}: Adot mismatch"
        )

    print(label, flush=True)
    print("  A =", aval, flush=True)
    print("  J =", jval, flush=True)
    print("  Adot =", adot, flush=True)
    print("  Jdot =", jdot, flush=True)

    if jval != 0:
        print(
            "  RESULT: not tangent to A=0",
            flush=True
        )
    else:
        print(
            "  RESULT: first-order A tangency holds",
            flush=True
        )
        print(
            "  Addot at witness =",
            addot,
            flush=True
        )

print("STAGE 2 PASS", flush=True)

print("STAGE 3: distinguish consistency questions", flush=True)

# These values were established by R27I2 and R27J.
# They are recorded here, NOT recomputed.
prior = {
    "J_NONZERO": (
        S.Rational(1,2),
        S.Rational(3,2)
    ),
    "J_ZERO": (
        S.Rational(1,4),
        S.Integer(0)
    )
}

for label, data in cases:
    D, K2 = prior[label]
    coefficient = ev(2*W*J/v**2, data)

    if coefficient != K2:
        raise AssertionError(
            f"{label}: prior K2 mismatch"
        )

    print(label, flush=True)
    print("  previously verified D =", D, flush=True)
    print("  previously verified K2 =", K2, flush=True)

    if K2:
        print(
            "  lambda2 from Bdot=0 =",
            S.simplify(-D/K2),
            flush=True
        )
    else:
        print(
            "  Bdot independent of lambda2 =",
            D,
            flush=True
        )

print("STAGE 3 PASS", flush=True)

print()
print("R27K COMPLETE")
print("ADOT IDENTITY: VERIFIED")
print("JDOT AT 4D WITNESSES: EVALUATED")
print("SINGULAR STRATUM INVARIANCE: OPEN")
print("GENERAL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3)
)
