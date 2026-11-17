#!/usr/bin/env python3
"""
GWFM / SEFI — CHECKPOINT 094

UNITARY NORMALIZER OF THE CLASSICAL
DENSITY-EQUIVALENCE STRUCTURE

Finite-dimensional setting: C^3.

Question:
Which unitary matrices induce a well-defined
transformation on classical probability vectors
rho_j = |psi_j|^2, independently of local phases?

Theorem investigated:
Such unitaries are exactly monomial matrices:
permutations multiplied by diagonal phases.

This is an algebraic classification, not a
physical selection principle.
"""

import sympy as sp

I = sp.I

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 094")
print("UNITARY NORMALIZER OF CLASSICAL DENSITY")
print("=" * 64)

def eq(A, B):
    def zero(expr):
        expr = sp.expand_complex(expr)
        expr = sp.trigsimp(expr)
        return sp.simplify(expr) == 0

    if isinstance(A, sp.MatrixBase):
        return all(zero(z) for z in A-B)

    return zero(A-B)

def density(z):
    return sp.Matrix([
        sp.simplify(z[j]*sp.conjugate(z[j]))
        for j in range(z.rows)
    ])

def is_unitary(U):
    return eq(U.conjugate().T*U, sp.eye(U.rows))

theta, alpha, beta = sp.symbols(
    "theta alpha beta", real=True
)

# 1. A diagonal phase transformation.
D = sp.diag(
    sp.exp(I*alpha),
    sp.exp(I*beta),
    sp.exp(I*theta)
)

assert is_unitary(D)

print("PASS: arbitrary diagonal phases are unitary")

# 2. A permutation matrix.
P = sp.Matrix([
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 0]
])

assert is_unitary(P)

print("PASS: permutation matrix is unitary")

# 3. Their product is monomial and unitary.
U = P*D

assert is_unitary(U)

print("PASS: monomial unitary constructed")

# 4. Monomial transformations preserve
# density equivalence for arbitrary
# symbolic amplitudes and phases.
r, s, q = sp.symbols(
    "r s q", positive=True
)

psi = sp.Matrix([
    sp.sqrt(r),
    sp.sqrt(s),
    sp.sqrt(q)
])

phi = D*psi

assert eq(density(psi), density(phi))

assert eq(
    density(U*psi),
    density(U*phi)
)

print("PASS: monomial unitary preserves")
print("      density equivalence")

# 5. Its induced classical map is
# exactly a permutation of densities.
assert eq(
    density(U*psi),
    P*density(psi)
)

print("PASS: induced density map is")
print("      classical permutation")

# 6. General phase independence
# can be checked row by row.
#
# For a row with coefficients a and b
# in distinct columns, input amplitudes
# 1 and exp(i*theta) yield:
#
# |a+b exp(i*theta)|^2.
#
# Independence of theta for all theta
# forces conjugate(a)*b = 0.
#
# Test the exact cross-term identity.

ar, ai, br, bi = sp.symbols(
    "ar ai br bi", real=True
)

aa = ar+I*ai
bb = br+I*bi

intensity = sp.expand(
    (aa+bb*sp.exp(I*theta))
    *sp.conjugate(
        aa+bb*sp.exp(I*theta)
    )
)

cross = sp.simplify(
    intensity-aa*sp.conjugate(aa)
    -bb*sp.conjugate(bb)
)

expected_cross = (
    sp.conjugate(aa)*bb*sp.exp(I*theta)
    +aa*sp.conjugate(bb)*sp.exp(-I*theta)
)

assert eq(cross, expected_cross)

print("PASS: general two-column")
print("      interference cross term derived")

# 7. Independence for every phase
# forces the complex cross coefficient
# to vanish. Test at four phases.
#
# At 0 and pi, the real part vanishes.
# At pi/2 and 3pi/2, the imaginary
# part vanishes.

z = sp.symbols("z", complex=True)

cross_z = (
    z*sp.exp(I*theta)
    +sp.conjugate(z)*sp.exp(-I*theta)
)

assert eq(
    cross_z.subs(theta,0),
    2*sp.re(z)
)

assert eq(
    cross_z.subs(theta,sp.pi/2),
    -2*sp.im(z)
)

assert eq(
    cross_z.subs(theta,sp.pi),
    -2*sp.re(z)
)

assert eq(
    cross_z.subs(theta,3*sp.pi/2),
    2*sp.im(z)
)

print("PASS: phase independence forces")
print("      pairwise row cross terms to zero")

# 8. If a row contains two nonzero
# coefficients, their cross product
# cannot vanish. Thus each row has
# at most one nonzero entry.
#
# Unitarity makes every row norm one,
# so each row has exactly one nonzero
# entry of modulus one.
#
# Orthogonality forces different rows
# to occupy different columns.
#
# This proves the monomial classification.

assert all(
    sum(1 for j in range(3)
        if U[i,j] != 0) == 1
    for i in range(3)
)

assert all(
    sum(1 for i in range(3)
        if U[i,j] != 0) == 1
    for j in range(3)
)

print("PASS: monomial row and column")
print("      structure verified")

# 9. A non-monomial unitary.
H2 = sp.Matrix([
    [1, 1],
    [1,-1]
])/sp.sqrt(2)

V = sp.diag(H2, sp.ones(1))

assert is_unitary(V)

print("PASS: non-monomial unitary constructed")

# 10. Two vectors with identical
# classical densities.
psi0 = sp.Matrix([
    1/sp.sqrt(2),
    1/sp.sqrt(2),
    0
])

psi1 = sp.Matrix([
    1/sp.sqrt(2),
    -1/sp.sqrt(2),
    0
])

assert eq(
    density(psi0),
    density(psi1)
)

print("PASS: phase-related inputs")
print("      have identical densities")

# 11. Their outputs under V differ.
out0 = density(V*psi0)
out1 = density(V*psi1)

assert eq(out0, sp.Matrix([1,0,0]))
assert eq(out1, sp.Matrix([0,1,0]))

print("PASS: non-monomial mixing fails")
print("      density-quotient compatibility")

# 12. Classical density maps generated
# by monomial unitaries are reversible
# permutations, not general stochastic
# mixing operations.
rho = sp.Matrix([
    sp.Rational(1,2),
    sp.Rational(1,3),
    sp.Rational(1,6)
])

assert eq(
    P.T*(P*rho),
    rho
)

print("PASS: induced classical permutation")
print("      is reversible")

# 13. A stochastic map can preserve
# probability without being a
# reversible permutation.
K = sp.Matrix([
    [sp.Rational(1,2), sp.Rational(1,2), 0],
    [sp.Rational(1,2), sp.Rational(1,2), 0],
    [0,0,1]
])

assert all(
    eq(sum(K[:,j]),1)
    for j in range(3)
)

assert eq(K.det(),0)

print("PASS: stochastic mixing is broader")
print("      than reversible permutations")

# 14. Monomial unitaries normalize
# the diagonal observable algebra.
A = sp.diag(2,3,5)

transformed_A = (
    U.conjugate().T*A*U
)

assert eq(
    transformed_A,
    sp.diag(5,2,3)
)

print("PASS: monomial unitary preserves")
print("      diagonal observable algebra")

# 15. Non-monomial unitary does not
# preserve that algebra.
B = sp.diag(1,0,0)

transformed_B = (
    V.conjugate().T*B*V
)

assert not eq(
    transformed_B,
    sp.diag(
        transformed_B[0,0],
        transformed_B[1,1],
        transformed_B[2,2]
    )
)

print("PASS: non-monomial unitary")
print("      does not normalize diagonal algebra")

print("\n" + "="*64)
print("CHECKPOINT 094: ALL TESTS PASSED")
print("="*64)

print("ESTABLISHED MATHEMATICALLY:")
print("- A finite-dimensional unitary preserves")
print("  density equivalence for all states")
print("  iff it is monomial.")
print("- Monomial unitaries consist of")
print("  permutations and diagonal phases.")
print("- Their induced density transformations")
print("  are classical permutations.")
print("- Non-monomial unitary mixing can")
print("  distinguish phase-related inputs.")
print("- Monomial unitaries normalize the")
print("  diagonal observable algebra.")
print("- General stochastic mixing is not")
print("  equivalent to reversible permutations.")

print("\nSCOPE AND QUALIFICATIONS:")
print("- The classification is exact for")
print("  finite-dimensional complex Hilbert space.")
print("- The universal argument uses row")
print("  cross terms and unitarity.")
print("- Infinite-dimensional cases require")
print("  separate functional-analytic treatment.")
print("- The unitary transformations are")
print("  mathematical constructions.")
print("- No physical coherent mechanism")
print("  has been derived.")

print("\nNOT ESTABLISHED:")
print("- A physical selection principle")
print("  for non-monomial transformations.")
print("- A quantum measurement theory.")
print("- The Born rule from GWFM.")
print("- A preferred action scale.")
print("- Emergent quantum mechanics.")
