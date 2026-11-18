"""CP098-T19: Exact two-parameter closure-family test.

Tests the T16 identity for arbitrary real t and n
within the specified regular constrained family.

Not a proof on the general constraint surface.
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

# Load T17's verified definitions without repeating its
# diagnostic output. Its assertions still execute.
with contextlib.redirect_stdout(io.StringIO()):
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
required = env["required"]

u = env["u"]
q = env["q"]
P = env["P"]
p1 = env["p1"]
I = env["I"]

t, n = sp.symbols("t n", real=True)

# Keep the original geometry and constrained transverse
# momenta, but vary two independent canonical momenta.
point[ps] = t
point[p0[0]] = 85 - t**2/2
point[p0[1]] = n

print("=== CP098-T19 ===", flush=True)

for name, expr in [
    ("phi_T", phiT),
    ("phi_N", phiN),
    ("chi_T", chiT),
    ("chi_N", chi),
    ("psi_T", psi),
]:
    result = sp.simplify(at(expr))
    print(name, "=", result, flush=True)
    assert result == 0, f"{name} not zero"

print("PASS: constraints hold identically", flush=True)

print("Computing canonical contributions...", flush=True)

part_u = canonical_pair(psi, chi, u, p1)
print("u,p1 =", part_u, flush=True)

part_q = canonical_pair(psi, chi, q, P)
print("q,P =", part_q, flush=True)

part_I = canonical_pair(psi, chi, [I], [ps])
print("internal =", part_I, flush=True)

direct = sp.factor(part_u + part_q + part_I)
target = sp.factor(at(required))
difference = sp.factor(direct - target)

print("Direct =", direct, flush=True)
print("Required =", target, flush=True)
print("Difference =", difference, flush=True)

assert difference == 0, (
    "Symbolic closure condition failed"
)

print("PASS: identity for arbitrary real t,n")
print("LIMIT: fixed geometry and other canonical data")
print("GENERAL DIRAC CLOSURE: NOT ESTABLISHED")
