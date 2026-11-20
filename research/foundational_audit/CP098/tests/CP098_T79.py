"""CP098-T79: Point-2 canonical E and full tangency audit.

Diagnostic only. Do not infer PASS from exit status.
"""

import contextlib
import io
import runpy
from pathlib import Path

import numpy as np
import sympy as sp

print("=== CP098-T79 ===", flush=True)

source = Path(__file__).with_name("CP098_T75.py")

with contextlib.redirect_stdout(io.StringIO()) as captured:
    ns = runpy.run_path(str(source))

if "CANONICAL E_DOT NUMERICAL CROSS-CHECK: PASS" not in captured.getvalue():
    raise RuntimeError("STOP: T75 prerequisite failed")

geometry = ns["geometry"]
internal = ns["internal"]
flow = ns["flow"]
canonical_E = ns["canonical_E"]
phiT = ns["phiT"]
phiN = ns["phiN"]
chiN = ns["chi"]
psiT = ns["psi"]

# Exact point 2 from T77.
v, k, a = 1.05, 0.43, -0.21
r, s = -0.14, 0.29
x, y = -0.18, 0.24
n, b, c = -0.16, 0.08, 0.17
I, pI, LT = -0.30, 0.42, -0.09

z = np.zeros(22)
z[0:4] = [v, 0, 0, 0]
z[4:8] = [a, v*v*k, 0, 0]
z[8:12] = [0, 0, r, s]
z[20:22] = [I, pI]

A, B, C = internal(z)
R = r*r+s*s
S = r*x+s*y

z[12:16] = [0, -2*A*k/v, x, y]
z[16:20] = [
    A*k*k+v**4*k*k*R/(4*C)-B,
    n, b, c
]

DA = float(ns["DAs"].subs({
    ns["ii"]: I,
    ns["pp"]: pI
}))

FN = (
    -n-6*A*a*k/v**2
    -v**3*k*S/(2*C)
    -v**2*k*a*R/(2*C)
    +2*k*DA
)
LN = -FN*v**3/(2*A)

def chiT(state):
    _, _, T, _ = geometry(state)
    return float(state[12:16] @ T)

constraints = {
    "phi_T": phiT,
    "phi_N": phiN,
    "chi_T": chiT,
    "chi_N": chiN,
    "psi_T": psiT,
}

print("CONSTRAINT VALUES", flush=True)
for name, function in constraints.items():
    print(
        f"{name}: {function(z):+.12e}",
        flush=True
    )

print("CONSTRAINT DERIVATIVES", flush=True)

for h in (1e-4, 1e-5):
    dz = flow(z, LT, LN, h/2)

    print(f"step={h:g}", flush=True)

    for name, function in constraints.items():
        derivative = (
            function(z+h*dz)
            -function(z-h*dz)
        )/(2*h)

        print(
            f"  d({name})/dt = "
            f"{derivative:+.12e}",
            flush=True
        )

E_surface = (
    -2*A*a*k/v**3
    +2*DA*k/v
    -n/v
    -a*k*v*R/(2*C)
    -k*v**2*S/(2*C)
)

print("CANONICAL E VS SURFACE E", flush=True)
print(
    f"E_surface={E_surface:+.15e}",
    flush=True
)

for h in (0.004, 0.002, 0.001, 0.0005):
    E_canonical = canonical_E(z, h)

    print(
        f"h={h:.4g} "
        f"E_canonical={E_canonical:+.15e} "
        f"difference={E_canonical-E_surface:+.9e}",
        flush=True
    )

print("DIRECTIONAL E DERIVATIVE", flush=True)

for h in (0.002, 0.001, 0.0005):
    dt = 4*h
    dz = flow(z, LT, LN, h/2)

    numerical = (
        canonical_E(z+dt*dz, h)
        -canonical_E(z-dt*dz, h)
    )/(2*dt)

    print(
        f"h={h:.4g} "
        f"E_dot_canonical={numerical:+.12e}",
        flush=True
    )

print("T79 DIAGNOSTIC COMPLETE")
print("T77: UNRESOLVED")
print("GENERAL SYMBOLIC CANONICAL PROOF: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
