"""CP098-T18: Two-parameter constrained-point checks.

Reuses the verified T17 definitions and localized
canonical differentiation. Tests multiple exact
points, not a general symbolic identity.
"""

import runpy
import sympy as sp
from pathlib import Path

source = Path(
    "research/foundational_audit/CP098/tests/CP098_T17.py"
)

assert source.is_file()
env = runpy.run_path(str(source))

point = env["point"]
p0 = env["p0"]
ps = env["ps"]

at = env["at"]
canonical_pair = env["canonical_pair"]

psi = env["psi"]
chi = env["chi"]
phiT = env["phiT"]
phiN = env["phiN"]
chiT = env["chiT"]
F_N = env["F_N"]
required = env["required"]

u = env["u"]
q = env["q"]
P = env["P"]
p1 = env["p1"]
I = env["I"]

print("\n=== CP098-T18 ===", flush=True)

cases = [
    (sp.Integer(0), sp.Integer(0)),
    (sp.Integer(2), sp.Integer(3)),
    (sp.Rational(-1, 2), sp.Integer(-4)),
]

for t, n in cases:
    point[ps] = t
    point[p0[0]] = 85 - t**2/2
    point[p0[1]] = n

    assert all(
        at(c) == 0
        for c in (phiT, phiN, chiT, chi, psi)
    ), "Constraint failure"

    direct = (
        canonical_pair(psi, chi, u, p1)
        + canonical_pair(psi, chi, q, P)
        + canonical_pair(psi, chi, [I], [ps])
    )

    difference = sp.simplify(
        direct - at(required)
    )

    print(
        f"t={t}, n={n}: "
        f"direct={direct}, "
        f"required={at(required)}, "
        f"difference={difference}",
        flush=True
    )

    assert difference == 0, (
        "T16 condition failed: investigate"
    )

print("PASS: all three constrained points")
print("GENERAL CLOSURE: NOT YET ESTABLISHED")
