"""CP098-T59: Relation between E/d and geometric evolution.

Cross-check independently established regular-surface
results from T41, T50, and T58.

All identities here are restricted to the regular
five-constraint surface.

No off-surface extension of E is assumed.
No gauge generator is constructed.
"""

import sympy as sp

v, k = sp.symbols("v k", positive=True)
A, C = sp.symbols("A C", real=True, nonzero=True)

a, n, r, s, x, y, DA = sp.symbols(
    "a n r s x y DA", real=True
)

LN = sp.symbols("lambda_N", real=True)

print("=== CP098-T59 ===", flush=True)

R = r*r+s*s
Q = r*x+s*y

# T41: regular-surface Poisson matrix entry.
E = -(
    4*A*C*a*k
    -4*C*DA*k*v**2
    +2*C*n*v**2
    +a*k*R*v**4
    +k*v**5*Q
)/(2*C*v**3)

# T40: normal second-class bracket.
d = 2*A/v**3

# T34/T50: normal forcing.
FN = (
    -n
    -6*A*a*k/v**2
    -v**3*k*Q/(2*C)
    -v**2*k*a*R/(2*C)
    +2*k*DA
)

# T50: multiplier determined by preservation.
LN_fixed = -FN/d

# T58: geometric evolution on regular surface,
# before the multiplier is fixed.
tdot = -2*a*k + LN/v

tdot_fixed = sp.factor(
    tdot.subs(LN, LN_fixed)
)

ratio = sp.factor(E/d)

print("E/d =", ratio, flush=True)
print("Fixed lambda_N =", sp.factor(LN_fixed),
      flush=True)
print("Fixed t_dot =", tdot_fixed, flush=True)

# First independent cross-check:
# E differs from FN/v by a simple geometric term.
identity_1 = sp.factor(
    E-FN/v-4*A*a*k/v**3
)

print("E - FN/v - 4*A*a*k/v^3 =",
      identity_1, flush=True)

assert identity_1 == 0
print("PASS: E/FN relation", flush=True)

# Second cross-check:
# E/d = 2*a*k - lambda_N_fixed/v.
identity_2 = sp.factor(
    ratio-(2*a*k-LN_fixed/v)
)

print("Multiplier relation residual =",
      identity_2, flush=True)

assert identity_2 == 0
print("PASS: multiplier relation", flush=True)

# Main result.
identity_3 = sp.factor(
    tdot_fixed+ratio
)

print("t_dot_fixed + E/d =", identity_3,
      flush=True)

assert identity_3 == 0
print("PASS: geometric evolution identity",
      flush=True)

# The identity is conditional on the normal
# multiplier consistency equation.
#
# Before fixing lambda_N:
#
# t_dot + E/d = (lambda_N - LN_fixed)/v.

unfixed_residual = sp.factor(
    tdot+ratio-(LN-LN_fixed)/v
)

print("Unfixed multiplier residual =",
      unfixed_residual, flush=True)

assert unfixed_residual == 0

print("PASS: multiplier dependence retained",
      flush=True)

# Re-express the corrected tertiary first-class
# combination using the verified on-surface
# coefficient identity:
#
# G3 = psi_T + t*chi_N - (E/d)*phi_N
#
#    = psi_T + t*chi_N + t_dot_fixed*phi_N
#
# This equality is between coefficient values
# evaluated on the regular constraint surface.

phiN, chiN, psiT = sp.symbols(
    "phi_N chi_N psi_T", real=True
)

t = v*k

G3_original = (
    psiT+t*chiN-ratio*phiN
)

G3_rewritten = (
    psiT+t*chiN+tdot_fixed*phiN
)

basis_residual = sp.factor(
    G3_original-G3_rewritten
)

assert basis_residual == 0

print("PASS: corrected tertiary basis rewritten",
      flush=True)

print("REGULAR-SURFACE IDENTITY:")
print("t_dot_fixed = -E/d")
print("G3 = psi_T + t*chi_N + t_dot_fixed*phi_N")

print("GEOMETRIC COEFFICIENT RELATION: VERIFIED")
print("OFF-SURFACE EXTENSION OF E/d: OPEN")
print("FULL TIME-DEPENDENT GENERATOR: OPEN")
print("CANONICAL ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
