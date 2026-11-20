"""CP098-T49: Exact tangential tertiary-chain audit.

Reconstructs T48's unrestricted-momentum result
and compares it against the exact constraints.

Algebraic cross-check; not independent Poisson
differentiation.
"""

import sympy as sp

v, k = sp.symbols("v k", positive=True)
A, B, C = sp.symbols("A B C", nonzero=True, real=True)

P0, P1, P2, P3 = sp.symbols(
    "P0 P1 P2 P3", real=True
)
p10, p11, p00 = sp.symbols(
    "p10 p11 p00", real=True
)

P_sq = -P0**2 + P1**2 + P2**2 + P3**2

phi_T = P0
phi_N = P1
chi_T = p10
chi_N = p11 + 2*A*k/v

psi_T = (
    p00 + B - A*k**2
    - v**4*k**2*P_sq/(4*C)
)

# T48 direct canonical result.
bracket = (
    3*A*k**2 - B + k*p11*v - p00
    + v**4*k**2*P_sq/(4*C)
)

candidate = -psi_T + v*k*chi_N
residual = sp.factor(bracket-candidate)

print("=== CP098-T49 ===", flush=True)
print("Exact identity residual =", residual, flush=True)

assert residual == 0
print("PASS: exact tangential chain identity", flush=True)

primary_secondary = {
    phi_T: 0,
    phi_N: 0,
    p10: 0,
    p11: -2*A*k/v
}

reduced = sp.factor(bracket.subs(primary_secondary))
expected = sp.factor((-psi_T).subs(primary_secondary))

print("Primary-secondary bracket =", reduced, flush=True)
print("Negative tertiary constraint =", expected, flush=True)

assert sp.factor(reduced-expected) == 0

# The tertiary condition is genuinely additional:
# its derivative with respect to p00 is nonzero.
assert sp.diff(psi_T, p00) == 1

print("PASS: tertiary constraint required", flush=True)

# No extra quadratic correction is needed.
incorrect_extra = v**4*k**2*(P0**2-P1**2)/(2*C)
assert sp.factor(bracket-(candidate+incorrect_extra)) != 0

print("PASS: redundant quadratic term rejected", flush=True)
print("TANGENTIAL TERTIARY CHAIN: VERIFIED")
print("NORMAL SECONDARY CHAIN: OPEN")
print("GAUGE GENERATOR: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
