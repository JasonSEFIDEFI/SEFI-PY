"""CP098-T34: Tangential multiplier consistency audit.

Checks the proposed tertiary identity and the
normal-multiplier contribution on the regular
external constrained family.

Internal-sector result is imported algebraically
from T33, not recomputed here.

No general closure claim.
"""

import sympy as sp

v, k = sp.symbols("v k", positive=True)
a, n, r, s, x, y = sp.symbols(
    "a n r s x y", real=True
)
A, C = sp.symbols("A C", real=True, nonzero=True)
DA = sp.symbols("DA", real=True)

R = r**2 + s**2
Q = r*x + s*y

# T29/T30 external F_N.
FN_ext = (
    -n
    - 6*A*a*k/v**2
    - v**3*k*Q/(2*C)
    - v**2*k*a*R/(2*C)
)

# T25 internal contribution.
FN = FN_ext + 2*k*DA

# T32 external bracket.
psi_chi_ext = (
    n/v
    + 2*A*a*k/v**3
    + a*k*v*R/(2*C)
    + k*v**2*Q/(2*C)
)

# T33 internal bracket.
psi_chi_int = -2*k*DA/v

psi_chi = sp.factor(
    psi_chi_ext + psi_chi_int
)

# T16 proposed identity.
Tq = -a

candidate = (
    4*A*k*Tq/v**3
    - FN/v
)

print("=== CP098-T34 ===", flush=True)

difference = sp.factor(psi_chi-candidate)

print("F_N =", sp.factor(FN))
print("{psi_T,chi_N} =", psi_chi)
print("Identity residual =", difference)

assert difference == 0

print("PASS: combined tertiary identity")

# Established normal bracket:
# {chi_N,phi_N}=+2A/v^3.
#
# Normal consistency:
# 0 = F_N + lambda_N*(2A/v^3).
#
# Hence lambda_N=-v^3*F_N/(2A).

lambda_N = -v**3*FN/(2*A)

normal_residual = sp.factor(
    FN + lambda_N*(2*A/v**3)
)

assert normal_residual == 0

print("PASS: normal multiplier equation")

# The complete tangential equation requires:
#
# {psi_T,H_c}
# + lambda_T*{psi_T,phi_T}
# + lambda_N*{psi_T,phi_N}.
#
# These brackets must be independently computed.
# Do not infer closure from the preceding identities.

print("NEXT REQUIRED BRACKETS:")
print("  {psi_T,H_c}")
print("  {psi_T,phi_T}")
print("  {psi_T,phi_N}")
print("TERTIARY CONSISTENCY: NOT YET RESOLVED")
print("DYNAMICAL CONSISTENCY: OPEN")
