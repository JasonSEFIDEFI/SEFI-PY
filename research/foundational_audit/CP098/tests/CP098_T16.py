"""CP098-T16: Tangent-tertiary closure dependency.

Establishes the exact algebraic condition that the
uncomputed general bracket {psi_T,chi_N} must satisfy.

This is a dependency analysis, not a closure proof.
"""

import sympy as sp

v, kappa, A = sp.symbols(
    "v kappa A", positive=True
)

F_N = sp.Symbol("F_N")
B_TN = sp.Symbol("B_TN")
Tq = sp.Symbol("Tq")

# Previously established tangent-secondary bracket.
B_TT = 4*A*kappa**2/v

# Constraint-surface Hamiltonian bracket.
H_bracket = (
    v**2*kappa*B_TN
    - Tq*B_TT
)

# Normal multiplier contribution after elimination.
multiplier_term = v*kappa*F_N

residual = sp.factor(
    H_bracket + multiplier_term
)

required_B_TN = sp.solve(
    residual, B_TN
)

assert len(required_B_TN) == 1

candidate = required_B_TN[0]

assert sp.simplify(
    residual.subs(B_TN, candidate)
) == 0

print("=== CP098-T16 ===")
print("Tangent-tertiary residual:")
print(residual)
print()
print("Required normal-secondary bracket:")
print("B_TN =", candidate)
print()
print("PASS: closure condition derived algebraically")
print("STATUS: general canonical bracket NOT verified")
print("NEXT: independently compute {psi_T,chi_N}")
