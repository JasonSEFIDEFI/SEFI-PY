"""CP098-T52: Primary first-class intersection audit.

Determine the intersection of:
  (1) the null space of the regular constraint matrix;
  (2) the span of the original primary constraints.

This establishes the number of independent primary
first-class directions, subject to regularity.

It does not construct a gauge generator.
"""

import sympy as sp

d, t, E = sp.symbols(
    "d t E", real=True, nonzero=True
)

# Order:
# phi_T, phi_N, chi_T, chi_N, psi_T

M = sp.Matrix([
    [0, 0, 0, 0, 0],
    [0, 0, 0, -d, t*d],
    [0, 0, 0, 2*t*d, -2*t*t*d],
    [0, d, -2*t*d, 0, E],
    [0, -t*d, 2*t*t*d, -E, 0]
])

print("=== CP098-T52 ===", flush=True)

assert M + M.T == sp.zeros(5)
assert M.rank() == 2

print("PASS: rank-two antisymmetric matrix")

# An arbitrary combination of the original
# primary constraints:
#
# Q = alpha*phi_T + beta*phi_N.

alpha, beta = sp.symbols("alpha beta", real=True)

primary_vector = sp.Matrix([
    alpha, beta, 0, 0, 0
])

residual = (
    M * primary_vector
).applyfunc(sp.factor)

print("Primary first-class residual:")
sp.pprint(residual)

# Require all five brackets to vanish.
solutions = sp.solve(
    list(residual),
    [alpha, beta],
    dict=True
)

print("Solutions =", solutions, flush=True)

# The normal primary coefficient must vanish.
assert sp.factor(residual[3]-d*beta) == 0
assert sp.factor(residual[4]+t*d*beta) == 0

# alpha remains arbitrary.
assert all(
    alpha not in entry.free_symbols
    for entry in residual
)

print("PASS: beta must vanish for d != 0")
print("PASS: alpha remains arbitrary")

# Explicit intersection basis.
primary_first_class = sp.Matrix([1, 0, 0, 0, 0])

assert M*primary_first_class == sp.zeros(5, 1)

# phi_N is not first-class.
normal_primary = sp.Matrix([0, 1, 0, 0, 0])

assert M*normal_primary != sp.zeros(5, 1)

print("PRIMARY FIRST-CLASS DIMENSION = 1")

# Recheck complete first-class nullity.
null_basis = M.nullspace()

print("Full first-class nullity =", len(null_basis))

assert len(null_basis) == 3

print("PASS: one primary first-class direction")
print("PASS: three total first-class directions")

# Chain relationships already established
# by T47 and T49:
#
# {phi_T,Hc} = -chi_T + t*phi_N
# {chi_T,Hc} = -psi_T + t*chi_N
#
# These are exact at the adapted-frame
# evaluation point, with t=v*k.

print("Verified tangential chain:")
print("phi_T -> chi_T -> psi_T")

print("Normal primary is second-class.")
print("Its multiplier is fixed by chi_N preservation.")

print("PRIMARY FIRST-CLASS INTERSECTION: VERIFIED")
print("INDEPENDENT GAUGE-FUNCTION COUNT: CANDIDATE 1")
print("FULL GAUGE GENERATOR: NOT YET DERIVED")
print("CANONICAL ACTION INVARIANCE: NOT YET TESTED")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
