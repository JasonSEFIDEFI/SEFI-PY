#!/usr/bin/env python3
"""
GWFM / SEFI — CHECKPOINT 090

CLASSIFICATION OF QUOTIENT-COMPATIBLE OBSERVABLES

Purpose:
Determine which operators preserve expectation
values under arbitrary local phase changes.

Exact result investigated:
In a finite-dimensional phase-space model,
an operator commuting with every diagonal
phase transformation must be diagonal.

Scope:
Finite-dimensional algebraic proof and
explicit counterexamples.

The corresponding continuum statement is
identified separately, not proved by the
finite-dimensional computation.

No quantum measurement assumptions.
"""

import sympy as sp

I = sp.I
n = 4

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 090")
print("CLASSIFICATION OF QUOTIENT-COMPATIBLE OBSERVABLES")
print("=" * 64)

def eq(a, b):
    if isinstance(a, sp.MatrixBase):
        return all(
            sp.simplify(v) == 0
            for v in a-b
        )
    return sp.simplify(a-b) == 0

# General complex matrix.
A = sp.Matrix(
    n, n,
    lambda i,j: sp.Symbol(
        f"A_{i}{j}", complex=True
    )
)

theta = sp.symbols(
    "theta0:4", real=True
)

U = sp.diag(
    *[sp.exp(I*t) for t in theta]
)

# 1. Local phases are unitary.
assert eq(
    U.conjugate().T*U,
    sp.eye(n)
)

print("PASS: arbitrary diagonal phase")
print("      transformations are unitary")

# 2. General matrix conjugation.
transformed = U.conjugate().T*A*U

for i in range(n):
    for j in range(n):
        expected = (
            sp.exp(I*(theta[j]-theta[i]))
            *A[i,j]
        )
        assert eq(transformed[i,j],expected)

print("PASS: exact phase transformation")
print("      of every matrix element")

# 3. Diagonal entries are invariant.
for i in range(n):
    assert eq(
        transformed[i,i],
        A[i,i]
    )

print("PASS: diagonal entries invariant")

# 4. Off-diagonal entries acquire
# independently variable phases.
for i in range(n):
    for j in range(n):
        if i == j:
            continue

        phases = {
            theta[i]: 0,
            theta[j]: sp.pi
        }

        for q in theta:
            if q not in phases:
                phases[q] = 0

        changed = sp.simplify(
            (transformed[i,j]-A[i,j])
            .subs(phases)
        )

        assert eq(changed,-2*A[i,j])

print("PASS: arbitrary local phases force")
print("      invariant off-diagonal entries")
print("      to vanish")

# 5. Every diagonal operator is invariant.
D = sp.diag(1,2,3,4)

assert eq(
    U.conjugate().T*D*U,
    D
)

print("PASS: diagonal multiplication")
print("      operators satisfy invariance")

# 6. The finite-dimensional commutant
# of all diagonal phase matrices is
# exactly the diagonal matrix algebra.
#
# Proof:
# - Every diagonal matrix commutes with U.
# - For i != j, a relative phase pi
#   changes A_ij to -A_ij.
# - Invariance requires A_ij = 0.
#
# This is a universal algebraic argument
# for the n-dimensional matrix model.

print("PASS: exact finite-dimensional")
print("      commutant classification")

# 7. The diagonal algebra is commutative.
D2 = sp.diag(5,7,11,13)

assert eq(
    D*D2-D2*D,
    sp.zeros(n)
)

print("PASS: quotient-compatible matrix")
print("      algebra is commutative")

# 8. Define a normalized vector and
# its rank-one projector.
psi = sp.Matrix([1,1,0,0])/sp.sqrt(2)

Q = psi*psi.conjugate().T

assert eq(
    (psi.conjugate().T*psi)[0],
    1
)

assert eq(Q*Q,Q)
assert eq(Q.conjugate().T,Q)

print("PASS: non-diagonal rank-one")
print("      orthogonal projector constructed")

# 9. Choose a relative phase reversal.
V = sp.diag(1,-1,1,1)

phi = V*psi

assert eq(
    V.conjugate().T*V,
    sp.eye(n)
)

# Classical discrete densities coincide.
for j in range(n):
    assert eq(
        sp.conjugate(psi[j])*psi[j],
        sp.conjugate(phi[j])*phi[j]
    )

print("PASS: distinct vectors have")
print("      identical discrete densities")

# 10. All diagonal expectations agree.
for observable in (D,D2,D+D2,D*D2):
    left = (
        psi.conjugate().T*observable*psi
    )[0]

    right = (
        phi.conjugate().T*observable*phi
    )[0]

    assert eq(left,right)

print("PASS: all tested diagonal")
print("      expectations agree")

# 11. Projector expectation differs.
q_psi = (
    psi.conjugate().T*Q*psi
)[0]

q_phi = (
    phi.conjugate().T*Q*phi
)[0]

assert eq(q_psi,1)
assert eq(q_phi,0)

print("PASS: non-diagonal projector")
print("      fails quotient compatibility")

# 12. State-only phase changes differ
# from simultaneous transformations.
Qprime = V*Q*V.conjugate().T

assert eq(
    (
        phi.conjugate().T*Qprime*phi
    )[0],
    q_psi
)

assert not eq(Qprime,Q)

print("PASS: consistent representation")
print("      change preserves expectation")

# 13. Expectation invariance for every
# vector implies operator invariance
# for self-adjoint operators.
#
# For B = U† A U - A self-adjoint:
# <z,Bz>=0 for every z implies B=0.
#
# In finite dimensions this follows
# from basis vectors and their sums.
#
# Verify the reconstruction identities
# on a general Hermitian 2x2 matrix.

r, s, h, j = sp.symbols(
    "r s h j", real=True
)

B = sp.Matrix([
    [r,h+I*j],
    [h-I*j,s]
])

e1 = sp.Matrix([1,0])
e2 = sp.Matrix([0,1])

def quadratic(z):
    return sp.expand(
        (z.conjugate().T*B*z)[0]
    )

assert eq(quadratic(e1),r)
assert eq(quadratic(e2),s)

assert eq(
    quadratic(e1+e2)-r-s,
    2*h
)

assert eq(
    quadratic(e1+I*e2)-r-s,
    -2*j
)

print("PASS: Hermitian operator determined")
print("      by its quadratic expectations")

# 14. Finite-dimensional conclusion.
assert eq(
    U.conjugate().T*D*U,
    D
)

assert not eq(
    V.conjugate().T*Q*V,
    Q
)

print("PASS: diagonal observables descend;")
print("      tested non-diagonal projector does not")

print("\n" + "=" * 64)
print("CHECKPOINT 090: ALL TESTS PASSED")
print("=" * 64)

print("ESTABLISHED MATHEMATICALLY:")
print("- In the finite-dimensional model,")
print("  local phases act diagonally.")
print("- An operator invariant under all")
print("  independent local phases is diagonal.")
print("- Every diagonal operator is invariant.")
print("- The resulting invariant algebra")
print("  is commutative.")
print("- A non-diagonal rank-one projector")
print("  distinguishes density-equivalent states.")
print("- Consistent unitary transformation of")
print("  state and operator preserves expectations.")

print("\nCONTINUUM CONNECTION:")
print("- On L2 of a standard phase space,")
print("  the bounded multiplication algebra")
print("  is a maximal abelian von Neumann algebra.")
print("- Its commutant consists of bounded")
print("  multiplication operators.")
print("- This continuum result is an established")
print("  functional-analytic theorem, not a")
print("  consequence proved by these finite tests.")
print("- Unbounded differential operators require")
print("  separate domain-sensitive treatment.")

print("\nNOT ESTABLISHED:")
print("- Physical selection of a larger algebra.")
print("- Quantum measurement or Born rule.")
print("- A preferred action scale.")
print("- Quantum dynamics emerging from GWFM.")
print("- A derivation of quantum mechanics.")
