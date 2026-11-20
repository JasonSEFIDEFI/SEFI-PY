"""CP098-T47: Primary-to-secondary chain identities.

Verify the exact adapted-frame off-surface identities
suggested by T46, using independent symbolic expressions.

This is an algebraic verification of T46, not an
independent Poisson differentiation.

No gauge-generator interpretation is assumed.
"""

import sympy as sp

v, k = sp.symbols("v k", positive=True)
A, C = sp.symbols("A C", real=True, nonzero=True)

P0, P1, P2, P3 = sp.symbols(
    "P0 P1 P2 P3", real=True
)

p10, p11 = sp.symbols("p10 p11", real=True)

print("=== CP098-T47 ===", flush=True)

# Adapted-frame constraints, without imposing them.
phi_T = P0
phi_N = P1
chi_T = p10
chi_N = p11 + 2*A*k/v

# Direct T46 expressions, retained off surface.
bracket_T = P1*k*v - p10

bracket_N = -(
    4*A*C*k
    - 2*C*P0*k*v**2
    + 2*C*p11*v
    + P0**2*k*v**4
    - P1**2*k*v**4
)/(2*C*v)

# Proposed exact chain identities.
chain_T = -chi_T + v*k*phi_N

chain_N = (
    -chi_N
    + v*k*phi_T
    - v**3*k*(phi_T**2-phi_N**2)/(2*C)
)

residual_T = sp.factor(bracket_T-chain_T)
residual_N = sp.factor(bracket_N-chain_N)

print("Tangential chain residual =", residual_T,
      flush=True)
print("Normal chain residual =", residual_N,
      flush=True)

assert residual_T == 0
assert residual_N == 0

print("PASS: exact primary-chain identities",
      flush=True)

# Restrict only to the primary constraint surface.
primary_surface = {P0: 0, P1: 0}

secondary_T = sp.factor(
    bracket_T.subs(primary_surface)
)
secondary_N = sp.factor(
    bracket_N.subs(primary_surface)
)

print("Primary-surface phi_T evolution =",
      secondary_T, flush=True)
print("Primary-surface phi_N evolution =",
      secondary_N, flush=True)

assert sp.factor(secondary_T+chi_T) == 0
assert sp.factor(secondary_N+chi_N) == 0

print("PASS: secondary constraints follow",
      flush=True)

# Check that the quadratic correction has no
# linear contribution on the primary surface.
quadratic = -v**3*k*(P0**2-P1**2)/(2*C)

assert quadratic.subs(primary_surface) == 0
assert sp.diff(quadratic, P0).subs(primary_surface) == 0
assert sp.diff(quadratic, P1).subs(primary_surface) == 0

print("PASS: quadratic correction identified",
      flush=True)

# The full off-surface identities must hold
# for arbitrary P2 and P3. They do not appear
# in the final expressions.
for z in (P2, P3):
    assert sp.diff(bracket_T, z) == 0
    assert sp.diff(bracket_N, z) == 0

print("PASS: transverse momenta unrestricted",
      flush=True)

print("PRIMARY-TO-SECONDARY CHAINS: VERIFIED")
print("SECONDARY-TO-TERTIARY CHAINS: OPEN")
print("GAUGE GENERATOR: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
