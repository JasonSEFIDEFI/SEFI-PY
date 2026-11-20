"""CP098-T76: Convergence-order audit of T75.

Reuses T75's canonical functions and regular point.
Checks second-order finite-difference convergence
and Richardson extrapolation.

Does not constitute a general symbolic proof.
"""

import contextlib
import io
import runpy
from pathlib import Path

import numpy as np

print("=== CP098-T76 ===", flush=True)

source = Path(__file__).with_name("CP098_T75.py")

if not source.is_file():
    raise SystemExit("STOP: T75 source missing")

# Execute the existing T75 checks before reusing
# its definitions. Preserve its source unchanged.
capture = io.StringIO()

with contextlib.redirect_stdout(capture):
    ns = runpy.run_path(str(source))

if "CANONICAL E_DOT NUMERICAL CROSS-CHECK: PASS" not in capture.getvalue():
    raise AssertionError("STOP: T75 prerequisite did not pass")

print("PASS: T75 prerequisite rerun", flush=True)

z = ns["z"]
LT = ns["LT"]
LN = ns["LN"]

flow = ns["flow"]
canonical_E = ns["canonical_E"]
expected = ns["expected"]

steps = [0.004, 0.002, 0.001, 0.0005]
estimates = []
errors = []

for h in steps:
    dt = 4*h
    dz = flow(z, LT, LN, h/2)

    Ep = canonical_E(z+dt*dz, h)
    Em = canonical_E(z-dt*dz, h)

    estimate = (Ep-Em)/(2*dt)
    error = abs(estimate-expected)

    estimates.append(estimate)
    errors.append(error)

    print(
        f"h={h:.6f} "
        f"estimate={estimate:.12f} "
        f"error={error:.6e}",
        flush=True
    )

# Require improvement at each halving.
for previous, current in zip(errors, errors[1:]):
    if not current < previous:
        raise AssertionError(
            "STOP: error did not decrease monotonically"
        )

print("PASS: monotonic error reduction", flush=True)

orders = []

for i in range(len(errors)-1):
    order = np.log2(errors[i]/errors[i+1])
    orders.append(order)

    print(
        f"observed_order_{i+1}={order:.6f}",
        flush=True
    )

# Central differences should approach order 2
# before floating-point cancellation dominates.
if not all(1.5 < p < 2.5 for p in orders):
    raise AssertionError(
        "STOP: observed convergence order inconsistent"
    )

print("PASS: second-order convergence", flush=True)

# Richardson extrapolation from the two
# smallest step sizes:
#
# D(h) = D(0) + c*h^2 + O(h^4)
#
# D_R = (4*D(h/2)-D(h))/3

richardson = (
    4*estimates[-1]-estimates[-2]
)/3

richardson_error = abs(richardson-expected)

print(
    f"Richardson_estimate={richardson:.12f}",
    flush=True
)
print(
    f"Richardson_error={richardson_error:.6e}",
    flush=True
)

if not richardson_error < errors[-1]:
    raise AssertionError(
        "STOP: extrapolation did not improve accuracy"
    )

print("PASS: Richardson improvement", flush=True)

print("NUMERICAL CONVERGENCE AUDIT: VERIFIED")
print("SECOND-ORDER BEHAVIOR: OBSERVED")
print("GENERAL SYMBOLIC CANONICAL PROOF: OPEN")
print("INDEPENDENT MULTI-POINT AUDIT: OPEN")
print("JACOBI CROSS-CHECK: OPEN")
print("GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
