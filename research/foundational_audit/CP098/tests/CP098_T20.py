"""CP098-T20: Geometry-dependent tangent closure checks.

Uses exact arithmetic, symbolic t,n, and multiple
regular adapted geometries.

No claim of general constraint-surface closure.
"""

import contextlib
import io
import runpy
from pathlib import Path
import sympy as sp

source = Path(
    "research/foundational_audit/CP098/tests/CP098_T17.py"
)
assert source.is_file()

with contextlib.redirect_stdout(io.StringIO()):
    env = runpy.run_path(str(source))

point = env["point"]
u, q = env["u"], env["q"]
P, p1, p0 = env["P"], env["p1"], env["p0"]
I, ps = env["I"], env["ps"]

at = env["at"]
bracket = env["canonical_pair"]

constraints = [
    env["phiT"],
    env["phiN"],
    env["chiT"],
    env["chi"],
    env["psi"],
]

t, n = sp.symbols("t n", real=True)

def set_vector(vec, values):
    for symbol, value in zip(vec, values):
        point[symbol] = sp.sympify(value)

def run_case(speed, curvature, longitudinal):
    v = sp.Rational(speed)
    k = sp.Rational(curvature)
    a = sp.Rational(longitudinal)

    assert v > 0 and k > 0

    # Timelike adapted frame:
    # T=(1,0,0,0), N=(0,1,0,0)
    # q0=a, q1=v^2*k
    set_vector(u, [v, 0, 0, 0])
    set_vector(q, [a, v**2*k, 0, 0])

    set_vector(P, [0, 0, 1, 0])
    set_vector(p1, [0, -2*k/v, 0, 0])

    # psi_T=0 with A=C=1, B=t^2/2.
    p00 = k**2 + v**4*k**2/4 - t**2/2
    set_vector(p0, [p00, n, 0, 0])

    point[I] = sp.Integer(1)
    point[ps] = t

    for name, expression in zip(
        ["phi_T", "phi_N", "chi_T", "chi_N", "psi_T"],
        constraints,
    ):
        value = sp.factor(at(expression))
        assert value == 0, (
            f"Constraint {name} failed: {value}"
        )

    psi = env["psi"]
    chi = env["chi"]

    part_u = bracket(psi, chi, u, p1)
    part_q = bracket(psi, chi, q, P)
    part_I = bracket(psi, chi, [I], [ps])

    direct = sp.factor(part_u + part_q + part_I)
    required = sp.factor(at(env["required"]))
    difference = sp.factor(direct - required)

    print(
        f"v={v}, kappa={k}, q0={a}",
        flush=True,
    )
    print("  direct =", direct, flush=True)
    print("  required =", required, flush=True)
    print("  difference =", difference, flush=True)

    assert difference == 0, (
        "Closure identity failed in this geometry"
    )

print("=== CP098-T20 ===", flush=True)

cases = [
    (2, 1, 0),
    (3, 2, 1),
    (4, 1, -2),
]

for case in cases:
    run_case(*case)

print("PASS: all geometry-dependent families")
print("Each result holds for arbitrary real t,n")
print("GENERAL DIRAC CLOSURE: NOT ESTABLISHED")
