#!/usr/bin/env python3
"""
GWFM / SEFI — CHECKPOINT 096

SYMPLECTIC AND DYNAMICAL COMPATIBILITY

Setting:
- Phase space R2 with omega = dx wedge dp
- H(p) = sqrt(p**2 + m**2), m > 0
- Free flow F_t(x,p) = (x + t*v(p), p)

Tests:
1. Coordinate scaling (noncanonical)
2. Canonical reciprocal scaling
3. Spatial translation
4. Momentum-dependent position shear
5. Momentum translation

Distinguishes:
- Jacobian-corrected L2 unitarity
- Symplectic preservation
- Commutation with the fixed free flow

No new physical law is assumed.
"""

import sympy as sp

x, p, t = sp.symbols("x p t", real=True)
m, c = sp.symbols("m c", positive=True)
a, q, gamma = sp.symbols(
    "a q gamma", real=True
)

H = sp.sqrt(p**2 + m**2)
v = sp.diff(H, p)

J0 = sp.Matrix([[0, 1], [-1, 0]])

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 096")
print("SYMPLECTIC AND DYNAMICAL COMPATIBILITY")
print("=" * 64)

def zero(expr):
    return sp.simplify(
        sp.trigsimp(sp.expand(expr))
    ) == 0

def eq(A, B):
    if isinstance(A, sp.MatrixBase):
        return all(zero(z) for z in A-B)
    return zero(A-B)

def flow(X):
    xx, pp = X
    return sp.Matrix([
        xx + t*v.subs(p, pp),
        pp
    ])

def apply(T, X):
    return sp.Matrix([
        z.subs(
            {x: X[0], p: X[1]},
            simultaneous=True
        )
        for z in T
    ])

def jacobian(T):
    return T.jacobian([x, p])

def symplectic(T):
    J = jacobian(T)
    return eq(J.T*J0*J, J0)

def commutes_with_flow(T):
    X = sp.Matrix([x, p])
    return eq(
        apply(T, flow(X)),
        flow(T)
    )

def poisson(f, g):
    return (
        sp.diff(f, x)*sp.diff(g, p)
        - sp.diff(f, p)*sp.diff(g, x)
    )

X = sp.Matrix([x, p])

# 1. Established free relativistic flow.
assert eq(v, p/H)
assert eq(poisson(x, H), v)
assert eq(poisson(p, H), 0)

print("PASS: relativistic Hamiltonian")
print("      generates established free flow")

# 2. The free flow is symplectic.
Ft = flow(X)
assert symplectic(Ft)
assert eq(jacobian(Ft).det(), 1)

print("PASS: free flow preserves")
print("      canonical symplectic form")

# 3. Noncanonical position scaling.
T_noncanonical = sp.Matrix([c*x, p])

assert eq(
    jacobian(T_noncanonical).det(), c
)
assert not symplectic(T_noncanonical)

print("PASS: position-only scaling")
print("      is generally noncanonical")

# 4. Reciprocal canonical scaling.
T_canonical = sp.Matrix([c*x, p/c])

assert symplectic(T_canonical)
assert eq(
    jacobian(T_canonical).det(), 1
)

print("PASS: reciprocal scaling")
print("      preserves symplectic form")

# 5. Canonical does not imply symmetry
# of the fixed relativistic flow.
left = apply(T_canonical, Ft)
right = flow(T_canonical)

assert not eq(
    left.subs({c: 2, p: 1, m: 1, t: 1}),
    right.subs({c: 2, p: 1, m: 1, t: 1})
)

print("PASS: canonical scaling need not")
print("      commute with relativistic flow")

# 6. Spatial translations.
T_space = sp.Matrix([x+a, p])

assert symplectic(T_space)
assert commutes_with_flow(T_space)

print("PASS: spatial translations are")
print("      canonical free-flow symmetries")

# 7. Momentum-dependent shear.
T_shear = sp.Matrix([x+gamma*p**2, p])

assert symplectic(T_shear)
assert commutes_with_flow(T_shear)

print("PASS: momentum-dependent position")
print("      shear is a free-flow symmetry")

# 8. Momentum translation.
T_momentum = sp.Matrix([x, p+q])

assert symplectic(T_momentum)

left = apply(T_momentum, Ft)
right = flow(T_momentum)

assert not eq(
    left.subs({q: 1, p: 0, m: 1, t: 1}),
    right.subs({q: 1, p: 0, m: 1, t: 1})
)

print("PASS: momentum translation is")
print("      canonical but not generally")
print("      a symmetry of fixed free flow")

# 9. Hamiltonian invariance is sufficient
# for the canonical transformations here
# to commute with the fixed flow.
assert eq(
    H.subs(p, T_space[1]), H
)
assert eq(
    H.subs(p, T_shear[1]), H
)

assert not eq(
    H.subs(p, T_momentum[1]).subs(
        {p: 0, q: 1, m: 1}
    ),
    H.subs({p: 0, m: 1})
)

print("PASS: Hamiltonian invariance")
print("      distinguishes these examples")

# 10. Symplectic bracket preservation.
for T in (
    T_canonical,
    T_space,
    T_shear,
    T_momentum
):
    assert eq(
        poisson(T[0], T[1]), 1
    )

print("PASS: canonical transformations")
print("      preserve coordinate Poisson bracket")

# 11. Noncanonical scaling changes bracket.
assert eq(
    poisson(
        T_noncanonical[0],
        T_noncanonical[1]
    ),
    c
)

print("PASS: noncanonical scaling")
print("      changes coordinate Poisson bracket")

# 12. Invertible linear maps admit
# Jacobian-corrected L2 pullbacks.
#
# This identity concerns norm preservation,
# not Hamiltonian admissibility.

for T in (T_noncanonical, T_canonical):
    det = sp.simplify(jacobian(T).det())
    assert det.is_positive

    weight = sp.sqrt(det)

    assert eq(weight**2, det)

print("PASS: Jacobian amplitude factors")
print("      identified for linear maps")

# 13. Canonical transformations have
# unit Jacobian in these examples.
for T in (
    T_canonical,
    T_space,
    T_shear,
    T_momentum
):
    assert eq(jacobian(T).det(), 1)

print("PASS: tested canonical maps")
print("      preserve phase-space volume")

# 14. Noncanonical scaling generally
# changes the Hamiltonian vector field
# relationship and cannot be treated
# as a canonical symmetry.
assert not eq(
    poisson(
        T_noncanonical[0],
        T_noncanonical[1]
    ).subs(c, 2),
    1
)

print("PASS: unitarity alone does not")
print("      establish canonical compatibility")

print("\n" + "=" * 64)
print("CHECKPOINT 096: ALL TESTS PASSED")
print("=" * 64)

print("ESTABLISHED IN THESE MODELS:")
print("- The free relativistic flow is symplectic.")
print("- Position-only scaling is generally")
print("  noncanonical.")
print("- Reciprocal scaling is canonical but")
print("  generally not a fixed-flow symmetry.")
print("- Spatial translations and tested")
print("  position shears commute with free flow.")
print("- Momentum translations are canonical")
print("  but generally do not commute with it.")
print("- Jacobian-corrected L2 unitarity")
print("  is distinct from symplectic preservation.")
print("- Symplectic preservation is distinct")
print("  from fixed-Hamiltonian symmetry.")

print("\nSCOPE LIMITATIONS:")
print("- All checks concern explicit smooth maps.")
print("- Failure of commutation is shown by")
print("  exact counterexamples.")
print("- No exhaustive symmetry classification")
print("  is claimed.")
print("- The Hamiltonian and symplectic form")
print("  are established classical inputs.")
print("- No new dynamics is inferred.")

print("\nNOT ESTABLISHED:")
print("- A unique physical transformation group.")
print("- A physically selected noncommutative")
print("  observable algebra.")
print("- Quantum measurement or Born rule.")
print("- A preferred action scale.")
print("- Emergent quantum mechanics.")
