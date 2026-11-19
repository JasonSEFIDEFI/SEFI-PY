import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28L ===", flush=True)

source = Path(__file__).with_name(
    "CP098_T90C_R28K2.py"
)
assert source.is_file()

# Reuse the verified witness and canonical bracket
# implementation, without changing R28K2.
with contextlib.redirect_stdout(io.StringIO()):
    previous = runpy.run_path(str(source))

ns = previous["ns"]
bracket_at = previous["bracket_at"]
at = previous["at"]
lambda2 = previous["lambda2"]

lam1 = S.symbols("lambda1", real=True)

assert S.simplify(
    lambda2 -
    7*S.sqrt(S.Integer(1187904394))
    /S.Integer(2878195)
) == 0

print("STAGE 1: witness regression", flush=True)

names = [
    "phi1", "phi2", "chi", "psi", "theta",
    "B", "A", "J", "Q", "pA"
]

for name in names:
    assert at(ns[name]) == 0, name

assert bracket_at(ns["B"], ns["H"]) == 0
print("PASS: ten conditions and D=0", flush=True)
print("lambda2 =", lambda2, flush=True)
print("STAGE 1 PASS", flush=True)

print("STAGE 2: simultaneous tangency", flush=True)

equations = []
unconditional_failures = []
lambda1_requirements = []

for name in names:
    f = ns[name]

    F = bracket_at(f, ns["H"])
    K1 = bracket_at(f, ns["phi1"])
    K2 = bracket_at(f, ns["phi2"])

    constant = S.simplify(F + lambda2*K2)
    residual = S.simplify(constant + lam1*K1)

    print()
    print("condition =", name, flush=True)
    print("F =", F, flush=True)
    print("K1 =", K1, flush=True)
    print("K2 =", K2, flush=True)
    print("dot =", residual, flush=True)

    equations.append(residual)

    if K1 == 0:
        if constant == 0:
            print("STATUS: satisfied", flush=True)
        else:
            print(
                "STATUS: OBSTRUCTION independent of lambda1",
                flush=True
            )
            unconditional_failures.append(name)
    else:
        required = S.simplify(-constant/K1)
        lambda1_requirements.append((name, required))
        print(
            "STATUS: requires lambda1 =",
            required,
            flush=True
        )

print()
print("STAGE 3: multiplier compatibility", flush=True)

if unconditional_failures:
    print(
        "UNCONDITIONAL OBSTRUCTIONS:",
        unconditional_failures,
        flush=True
    )
else:
    print(
        "PASS: no lambda1-independent obstruction",
        flush=True
    )

if lambda1_requirements:
    reference = lambda1_requirements[0][1]
    compatible = all(
        S.simplify(value-reference) == 0
        for _, value in lambda1_requirements
    )
    for name, value in lambda1_requirements:
        print(name, "requires", value, flush=True)
    print(
        "lambda1 requirements compatible =",
        compatible,
        flush=True
    )
else:
    compatible = True
    print(
        "lambda1 remains unrestricted by these equations",
        flush=True
    )

all_satisfied = (
    not unconditional_failures and compatible
)

print()
print("R28L COMPLETE", flush=True)
print(
    "SIMULTANEOUS TANGENCY COMPATIBLE =",
    all_satisfied,
    flush=True
)
print(
    "Ddot PRESERVATION: NOT YET CHECKED",
    flush=True
)
print(
    "FULL DIRAC CLOSURE: OPEN",
    flush=True
)
print(
    "PHYSICAL VALIDATION: NOT CLAIMED",
    flush=True
)
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3),
    flush=True
)
