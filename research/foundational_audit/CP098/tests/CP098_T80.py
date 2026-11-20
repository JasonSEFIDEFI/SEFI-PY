"""CP098-T80: Compare two off-surface extensions of E.

At the failing regular point, differentiate:
1. The T74 surface formula, evaluated using geometric scalars.
2. The original canonical bracket {chi_N,psi_T}.

Compare both with T74's analytical prediction.
Diagnostic only.
"""

import contextlib
import io
import runpy
from pathlib import Path
import numpy as np

print("=== CP098-T80 ===", flush=True)

source = Path(__file__).with_name("CP098_T75.py")

with contextlib.redirect_stdout(io.StringIO()) as captured:
    ns = runpy.run_path(str(source))

if "CANONICAL E_DOT NUMERICAL CROSS-CHECK: PASS" not in captured.getvalue():
    raise RuntimeError("STOP: T75 prerequisite failed")

geometry = ns["geometry"]
metric = ns["metric"]
internal = ns["internal"]
flow = ns["flow"]
canonical_E = ns["canonical_E"]

# Exact T77 point 2.
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

def internal_values(state):
    ii, pp = state[20:22]
    AA, BB, CC = internal(state)
    DA = (0.2+0.14*ii)*pp/1.4
    return AA, BB, CC, DA

DA = internal_values(z)[3]

FN = (
    -n-6*A*a*k/v**2
    -v**3*k*S/(2*C)
    -v**2*k*a*R/(2*C)
    +2*k*DA
)

LN = -FN*v**3/(2*A)

def E_surface_extension(state):
    vv, kk, T, N = geometry(state)
    q = state[4:8]
    P = state[8:12]
    p1 = state[12:16]
    p0 = state[16:20]

    aa = -metric(T, q)
    nn = float(p0 @ N)
    RR = metric(P, P)
    SS = metric(P, p1)

    AA, BB, CC, DD = internal_values(state)

    return (
        -2*AA*aa*kk/vv**3
        +2*DD*kk/vv
        -nn/vv
        -aa*kk*vv*RR/(2*CC)
        -kk*vv**2*SS/(2*CC)
    )

# T74 analytical prediction from the T75
# reference expression, using point-2 values.
subs = {
    ns["ii"]: I,
    ns["pp"]: pI
}

DC = float(ns["DCs"].subs(subs))
H = float(ns["Hs"].subs(subs))
J = float(ns["Js"].subs(subs))
K = float(ns["Ks"].subs(subs))

values = {
    ns["vv"]: v,
    ns["kk"]: k,
    ns["aa"]: a,
    ns["nn"]: n,
    ns["RR"]: R,
    ns["SS"]: S,
    ns["AA"]: A,
    ns["BB"]: B,
    ns["CC"]: C,
    ns["DD"]: DA,
    ns["DDC"]: DC,
    ns["HH"]: H,
    ns["JJ"]: J,
    ns["KK"]: K,
    ns["ltt"]: LT,
    ns["lnn"]: LN
}

expected = float(ns["Edot_symbolic"].subs(values))

print(f"T74 prediction: {expected:.12f}", flush=True)
print(f"lambda_N: {LN:.12f}", flush=True)

print(
    "E_surface_extension at point:",
    f"{E_surface_extension(z):.12f}",
    flush=True
)

print("DIRECTIONAL DERIVATIVES", flush=True)

for h in (0.004, 0.002, 0.001, 0.0005):
    dt = 4*h
    dz = flow(z, LT, LN, h/2)

    zp = z+dt*dz
    zm = z-dt*dz

    surface_dot = (
        E_surface_extension(zp)
        -E_surface_extension(zm)
    )/(2*dt)

    canonical_dot = (
        canonical_E(zp, h)
        -canonical_E(zm, h)
    )/(2*dt)

    print(
        f"h={h:.4g} "
        f"surface_dot={surface_dot:.12f} "
        f"canonical_dot={canonical_dot:.12f} "
        f"surface_minus_T74={surface_dot-expected:+.6e} "
        f"canonical_minus_T74={canonical_dot-expected:+.6e}",
        flush=True
    )

print("T80 DIAGNOSTIC COMPLETE", flush=True)
print("T77: UNRESOLVED", flush=True)
print("GENERAL CANONICAL PROOF: OPEN", flush=True)
print("PHYSICAL VALIDATION: NOT ESTABLISHED", flush=True)
