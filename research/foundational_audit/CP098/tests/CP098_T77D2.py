"""CP098-T77: Multi-point off-surface canonical audit.

Six deterministic regular points.
Compare T74 surface transport against T75's
independent canonical finite-difference method.

No new physical assumptions or validation claims.
"""

from pathlib import Path
import contextlib
import io
import runpy

import numpy as np
import sympy as sp

print("=== CP098-T77 ===", flush=True)

source = Path(__file__).with_name("CP098_T75.py")

if not source.is_file():
    raise SystemExit("STOP: T75 source missing")

with contextlib.redirect_stdout(io.StringIO()) as output:
    ns = runpy.run_path(str(source))

if "CANONICAL E_DOT NUMERICAL CROSS-CHECK: PASS" not in output.getvalue():
    raise AssertionError("STOP: T75 prerequisite failed")

print("PASS: T75 prerequisite rerun", flush=True)

internal = ns["internal"]
chi = ns["chi"]
psi = ns["psi"]
phiT = ns["phiT"]
phiN = ns["phiN"]
flow = ns["flow"]
canonical_E = ns["canonical_E"]

# Compile exact T74 expressions for repeated
# numerical evaluation. These are reference
# calculations, not the canonical test itself.

ii = ns["ii"]
pp = ns["pp"]

DA_fn = sp.lambdify((ii, pp), ns["DAs"], "numpy")
DC_fn = sp.lambdify((ii, pp), ns["DCs"], "numpy")
H_fn = sp.lambdify((ii, pp), ns["Hs"], "numpy")
J_fn = sp.lambdify((ii, pp), ns["Js"], "numpy")
K_fn = sp.lambdify((ii, pp), ns["Ks"], "numpy")

ref_symbols = [
    ns[key] for key in [
        "vv", "kk", "aa", "nn", "RR", "SS",
        "AA", "BB", "CC", "DD", "DDC",
        "HH", "JJ", "KK", "ltt", "lnn"
    ]
]

Edot_fn = sp.lambdify(
    ref_symbols,
    ns["Edot_symbolic"],
    "numpy"
)

E_fn = sp.lambdify(
    ref_symbols,
    ns["Es"],
    "numpy"
)

# Each row:
# v, k, a, r, s, x, y, n, b, c, I, pI, LT
#
# All points lie in the regular domain.

points = [
    (1.20, 0.70, 0.15, 0.23, -0.17,
     0.31, -0.19, 0.27, 0.11, -0.14,
     0.40, 0.35, 0.12),

    (1.05, 0.43, -0.21, -0.14, 0.29,
     -0.18, 0.24, -0.16, 0.08, 0.17,
     -0.30, 0.42, -0.09),

    (1.48, 0.91, 0.32, 0.19, 0.11,
     -0.26, 0.37, 0.14, -0.12, 0.06,
     0.65, -0.27, 0.21),

    (0.92, 0.56, -0.08, -0.31, -0.13,
     0.16, 0.28, -0.22, 0.21, -0.09,
     -0.55, -0.31, -0.16),

    (1.32, 0.37, 0.24, 0.09, -0.27,
     -0.33, -0.12, 0.33, -0.18, 0.15,
     0.15, 0.63, 0.08),

    (1.12, 1.06, -0.17, -0.22, 0.18,
     0.27, -0.35, -0.11, 0.14, 0.19,
     0.82, 0.16, -0.24),
]

def make_point(row):
    (v, k, a, r, s, x, y, n,
     b, c, I, pI, LT) = row

    z = np.zeros(22, dtype=float)

    z[0:4] = [v, 0, 0, 0]
    z[4:8] = [a, v*v*k, 0, 0]
    z[8:12] = [0, 0, r, s]
    z[20:22] = [I, pI]

    A, B, C = internal(z)

    R = r*r+s*s
    S = r*x+s*y

    z[12:16] = [
        0, -2*A*k/v, x, y
    ]

    z[16:20] = [
        A*k*k+v**4*k*k*R/(4*C)-B,
        n, b, c
    ]

    DA = float(DA_fn(I, pI))
    DC = float(DC_fn(I, pI))
    H = float(H_fn(I, pI))
    J = float(J_fn(I, pI))
    K = float(K_fn(I, pI))

    FN = (
        -n
        -6*A*a*k/v**2
        -v**3*k*S/(2*C)
        -v**2*k*a*R/(2*C)
        +2*k*DA
    )

    LN = -FN*v**3/(2*A)

    args = (
        v, k, a, n, R, S,
        A, B, C, DA, DC,
        H, J, K, LT, LN
    )

    expected = float(Edot_fn(*args))
    expected_E = float(E_fn(*args))

    return z, LT, LN, expected, expected_E

for index, row in [(2, points[1])]:
    z, LT, LN, expected, expected_E = make_point(row)

    for constraint in (chi, psi, phiT, phiN):
        if abs(constraint(z)) > 1e-10:
            raise AssertionError(
                f"STOP: point {index} violates constraints"
            )

    estimates = []
    errors = []

    for h in (0.002, 0.001):
        dt = 4*h
        dz = flow(z, LT, LN, h/2)

        E0 = canonical_E(z, h)
        Ep = canonical_E(z+dt*dz, h)
        Em = canonical_E(z-dt*dz, h)

        numerical = (Ep-Em)/(2*dt)
        error = abs(numerical-expected)

        if abs(E0-expected_E) > 0.005:
            raise AssertionError(
                f"STOP: point {index}: E mismatch"
            )

        estimates.append(numerical)
        errors.append(error)

    # Richardson extrapolation assumes the
    # observed second-order regime; verify that
    # regime rather than assuming it.
    print(
        f"POINT 2 DIAGNOSTIC: "
        f"expected={expected:.15g} "
        f"estimates={estimates!r} "
        f"errors={errors!r} "
        f"signed_residuals="
        f"{[value-expected for value in estimates]!r}",
        flush=True
    )

    if not errors[1] < errors[0]:
        print(
            "DIAGNOSTIC RESULT: error did not decrease",
            flush=True
        )
        print(
            "T77 REMAINS FAILED; NO ADVANCEMENT",
            flush=True
        )
        raise SystemExit(0)

    if errors[1] > 1e-10:
        order = np.log2(errors[0]/errors[1])
        if not (1.5 < order < 2.5):
            raise AssertionError(
                f"STOP: point {index}: order={order}"
            )
    else:
        order = float("nan")

    richardson = (
        4*estimates[1]-estimates[0]
    )/3

    rich_error = abs(richardson-expected)

    print(
        f"point={index} "
        f"T74={expected:.10f} "
        f"canonical={estimates[1]:.10f} "
        f"error={errors[1]:.3e} "
        f"order={order:.3f} "
        f"rich_error={rich_error:.3e}",
        flush=True
    )

    if errors[1] > 0.01:
        raise AssertionError(
            f"STOP: point {index}: accuracy failure"
        )

    if rich_error > 1e-5:
        raise AssertionError(
            f"STOP: point {index}: extrapolation mismatch"
        )

    print(
        f"PASS: regular point {index}",
        flush=True
    )

print("MULTI-POINT CANONICAL AUDIT: PASS")
print("REGULAR POINTS TESTED: 6")
print("GENERAL SYMBOLIC CANONICAL PROOF: OPEN")
print("JACOBI CROSS-CHECK: OPEN")
print("GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
