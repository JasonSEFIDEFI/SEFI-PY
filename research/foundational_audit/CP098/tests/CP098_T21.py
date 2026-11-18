"""CP098-T21: Symbolic five-parameter consistency family.

Parameters:
v > 0, kappa > 0, a, t, n real.

Uses exact canonical differentiation from T17.
No general Dirac closure claim.
"""

import contextlib
import io
import runpy
import sympy as sp
from pathlib import Path

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

v, k, a = sp.symbols(
    "v k a", positive=True, finite=True
)
t, n = sp.symbols("t n", real=True)

def assign(vec, values):
    point.update(zip(vec, values))

# Adapted timelike geometry.
assign(u, [v, 0, 0, 0])
assign(q, [a, v**2*k, 0, 0])

# Transverse canonical momentum.
assign(P, [0, 0, 1, 0])

# Enforce chi_N=0.
assign(p1, [0, -2*k/v, 0, 0])

# Enforce psi_T=0.
p00 = k**2 + v**4*k**2/4 - t**2/2
assign(p0, [p00, n, 0, 0])

point[I] = sp.Integer(1)
point[ps] = t

print("=== CP098-T21 ===", flush=True)
print("Checking constraints...", flush=True)

constraints = [
    ("phi_T", env["phiT"]),
    ("phi_N", env["phiN"]),
    ("chi_T", env["chiT"]),
    ("chi_N", env["chi"]),
    ("psi_T", env["psi"]),
]

for name, expr in constraints:
    result = sp.simplify(at(expr))
    print(name, "=", result, flush=True)
    assert result == 0, f"{name} failed"

print("PASS: constraints", flush=True)

psi = env["psi"]
chi = env["chi"]

print("Calculating (u,p1) bracket...", flush=True)
part_u = sp.factor(bracket(psi, chi, u, p1))
print("u,p1 =", part_u, flush=True)

print("Calculating (q,P) bracket...", flush=True)
part_q = sp.factor(bracket(psi, chi, q, P))
print("q,P =", part_q, flush=True)

print("Calculating internal bracket...", flush=True)
part_I = sp.factor(bracket(psi, chi, [I], [ps]))
print("internal =", part_I, flush=True)

direct = sp.factor(part_u + part_q + part_I)
required = sp.factor(at(env["required"]))

difference = sp.factor(direct - required)

print("Direct =", direct, flush=True)
print("Required =", required, flush=True)
print("Difference =", difference, flush=True)

assert difference == 0

print("PASS: symbolic five-parameter family")
print("GENERAL CONSTRAINT-SURFACE CLOSURE: UNRESOLVED")
