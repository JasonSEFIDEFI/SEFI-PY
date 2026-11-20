"""CP098-T82C: Efficient exact Hamiltonian identity audit.

Verifies the CP098 Hamiltonian decomposition
using exact algebraic identities rather than
expanding the complete geometric expressions.

Does not establish the constraint Jacobi identity.
"""

import sympy as sp

print("=== CP098-T82C ===", flush=True)

v, k, A, B, C = sp.symbols(
    "v k A B C", real=True, nonzero=True
)

Tq, p0u, p1q, p1T, p1N, P2 = sp.symbols(
    "Tq p0u p1q p1T p1N P2",
    real=True
)

# Geometric momentum decomposition:
# p1q = v^2*k*p1N - Tq*p1T.
geometric_relation = (
    p1q - v**2*k*p1N + Tq*p1T
)

Hc = (
    p0u+p1q
    +v*(B+A*k**2)
    -v**5*k**2*P2/(4*C)
)

# p0.T = p0.u/v.
psiT = (
    p0u/v+B-A*k**2
    -v**4*k**2*P2/(4*C)
)

chiN = p1N+2*A*k/v
chiT = p1T

decomposition = (
    v*psiT
    +v**2*k*chiN
    -Tq*chiT
)

residual = sp.factor(
    Hc-decomposition
)

if sp.simplify(residual-geometric_relation) != 0:
    raise AssertionError(
        "STOP: unexpected Hamiltonian residual"
    )

print("PASS: Hamiltonian residual equals geometric relation")

# Impose only the previously proved
# geometric identity, not the constraints.
reduced = sp.simplify(
    residual.subs(
        p1q,
        v**2*k*p1N-Tq*p1T
    )
)

if reduced != 0:
    raise AssertionError(
        f"STOP: nonzero residual {reduced}"
    )

print("PASS: exact Hamiltonian decomposition")
print("PASS: no constraint surface imposed")
print("PASS: no numerical approximation")
print("ACTUAL CONSTRAINT JACOBI IDENTITY: OPEN")
print("GENERAL E_DOT CANONICAL PROOF: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T82C COMPLETE")
