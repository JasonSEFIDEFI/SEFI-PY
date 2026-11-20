"""CP098-T78: isolate the point-2 transport discrepancy.

Compare canonical directional derivatives against
the scalar transport equations used in T74.

Diagnostic only. No mathematical PASS is asserted.
"""

import contextlib
import io
import runpy
from pathlib import Path

import numpy as np

print("=== CP098-T78 ===", flush=True)

source = Path(__file__).with_name("CP098_T75.py")

with contextlib.redirect_stdout(io.StringIO()) as captured:
    ns = runpy.run_path(str(source))

if "CANONICAL E_DOT NUMERICAL CROSS-CHECK: PASS" not in captured.getvalue():
    raise RuntimeError("STOP: T75 prerequisite failed")

print("PASS: T75 prerequisite", flush=True)

geometry = ns["geometry"]
internal = ns["internal"]
metric = ns["metric"]
flow = ns["flow"]

# Exactly the point-2 parameters from T77.
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
DC = float(ns["DCs"].subs({
    ns["ii"]: I,
    ns["pp"]: pI
}))
H = float(ns["Hs"].subs({
    ns["ii"]: I,
    ns["pp"]: pI
}))
J = float(ns["Js"].subs({
    ns["ii"]: I,
    ns["pp"]: pI
}))
K = float(ns["Ks"].subs({
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

def scalars(state):
    vv, kk, T, N = geometry(state)
    AA, BB, CC = internal(state)

    qq = state[4:8]
    PP = state[8:12]
    p1 = state[12:16]
    p0 = state[16:20]

    ii, pp = state[20:22]
    dai = (0.2+0.14*ii)*pp/1.4

    return np.array([
        vv,
        kk,
        -metric(T, qq),
        float(p0 @ N),
        metric(PP, PP),
        metric(PP, p1),
        AA,
        CC,
        dai
    ])

names = ["v", "k", "a", "n", "R", "S",
         "A", "C", "DA"]

expected = np.array([
    a,
    LN/v**2-3*a*k/v,
    LT-v**3*k*k,
    v*k*z[16]-v**3*k*(b*r+c*s)/(2*C),
    -2*S-2*LN*R/(v**2*k),
    (
        -x*x-y*y-r*b-s*c
        -LN*S/(v**2*k)
        -LT*R/v
        +a*LN*R/(k*v**3)
        -A*k*k*v*v*R/C
    ),
    v*DA,
    v*DC,
    v*H+v*k*k*J-v**5*k*k*R*K/4
])

print(f"lambda_N={LN:.12g}", flush=True)
print(
    "scalar     predicted          canonical"
    "          difference",
    flush=True
)

for h in (1e-4, 1e-5):
    dz = flow(z, LT, LN, h/2)

    actual = (
        scalars(z+h*dz)-scalars(z-h*dz)
    )/(2*h)

    print(f"--- step={h:g} ---", flush=True)

    for name, pred, got in zip(
        names, expected, actual
    ):
        print(
            f"{name:4s} "
            f"{pred: .11e} "
            f"{got: .11e} "
            f"{got-pred: .4e}",
            flush=True
        )

    print(
        "maximum absolute discrepancy:",
        np.max(np.abs(actual-expected)),
        flush=True
    )

print("DIAGNOSTIC COMPLETE — T77 STILL FAILED")
print("GENERAL CANONICAL PROOF: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
