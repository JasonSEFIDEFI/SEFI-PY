#!/usr/bin/env python3
"""
GWFM / SEFI — CHECKPOINT 097

GLOBAL CANONICAL SYMMETRIES OF
FREE RELATIVISTIC MOTION

Phase space: R2
Symplectic form: dx wedge dp
Hamiltonian: H = sqrt(p**2 + m**2), m > 0

Theorem under examination:

Every smooth global symplectomorphism
commuting with the complete free flow
has the form

T(x,p) = (eps*x + g(p), eps*p),

where eps = +1 or -1 and g is smooth.

The analytic proof uses:
1. Symplectic flow conjugacy.
2. Hamiltonian equality up to a constant.
3. Global surjectivity and the range [m,infinity).
4. Smoothness and connectedness.
5. The symplectic Jacobian condition.

Symbolic tests check the resulting formulas.
"""

import sympy as sp

x, p, t, m = sp.symbols(
    "x p t m", real=True
)
m = sp.symbols("m", positive=True)

a, b, c = sp.symbols(
    "a b c", real=True
)

H = sp.sqrt(p**2 + m**2)
v = p/H

J0 = sp.Matrix([[0,1],[-1,0]])

print("="*64)
print("GWFM / SEFI — CHECKPOINT 097")
print("GLOBAL FREE RELATIVISTIC FLOW SYMMETRIES")
print("="*64)

def eq(A,B):
    if isinstance(A,sp.MatrixBase):
        return all(
            sp.simplify(z) == 0
            for z in A-B
        )
    return sp.simplify(A-B) == 0

def flow(X):
    return sp.Matrix([
        X[0]+t*v.subs(p,X[1]),
        X[1]
    ])

def apply(T,X):
    return sp.Matrix([
        z.subs(
            {x:X[0],p:X[1]},
            simultaneous=True
        )
        for z in T
    ])

def is_symplectic(T):
    J = T.jacobian([x,p])
    return eq(J.T*J0*J,J0)

def commutes(T):
    base = sp.Matrix([x,p])
    return eq(
        apply(T,flow(base)),
        flow(T)
    )

# 1. Hamiltonian and velocity.
assert eq(sp.diff(H,p),v)
assert eq(sp.diff(v,p),
          m**2/(p**2+m**2)**sp.Rational(3,2))

print("PASS: relativistic velocity")
print("      and derivative verified")

# 2. Velocity is odd and H is even.
assert eq(v.subs(p,-p),-v)
assert eq(H.subs(p,-p),H)

print("PASS: parity properties verified")

# 3. General smooth shear family.
g = sp.Function("g")(p)

T_plus = sp.Matrix([
    x+g,p
])

assert is_symplectic(T_plus)
assert commutes(T_plus)

print("PASS: arbitrary smooth")
print("      momentum-dependent shear")
print("      is a canonical flow symmetry")

# 4. Momentum reversal combined
# with position reversal.
T_minus = sp.Matrix([
    -x+g,-p
])

assert is_symplectic(T_minus)
assert commutes(T_minus)

print("PASS: simultaneous position")
print("      and momentum reversal")
print("      is a canonical flow symmetry")

# 5. Both branches preserve H.
for T in (T_plus,T_minus):
    assert eq(
        H.subs(p,T[1]),H
    )

print("PASS: both symmetry branches")
print("      preserve the Hamiltonian")

# 6. Both branches are invertible.
#
# For T_plus:
# inverse (X,P) -> (X-g(P),P)
#
# For T_minus:
# inverse (X,P) -> (-X+g(-P),-P)

X, P = sp.symbols("X P",real=True)

plus_inverse = sp.Matrix([
    X-g.subs(p,P),P
])

minus_inverse = sp.Matrix([
    -X+g.subs(p,-P),-P
])

assert eq(
    plus_inverse.subs(
        {X:T_plus[0],P:T_plus[1]},
        simultaneous=True
    ),
    sp.Matrix([x,p])
)

assert eq(
    minus_inverse.subs(
        {X:T_minus[0],P:T_minus[1]},
        simultaneous=True
    ),
    sp.Matrix([x,p])
)

print("PASS: explicit global")
print("      inverses verified")

# 7. General canonical transformation
# with P = p and X = X(x,p).
#
# The symplectic condition gives
# partial_x X = 1.

F = sp.Function("F")(x,p)

J_plus = sp.Matrix([
    [sp.diff(F,x),sp.diff(F,p)],
    [0,1]
])

assert eq(
    J_plus.T*J0*J_plus,
    sp.diff(F,x)*J0
)

print("PASS: P=p forces")
print("      partial_x X=1")

# 8. For P=-p, symplecticity
# gives partial_x X=-1.

J_minus = sp.Matrix([
    [sp.diff(F,x),sp.diff(F,p)],
    [0,-1]
])

assert eq(
    J_minus.T*J0*J_minus,
    -sp.diff(F,x)*J0
)

print("PASS: P=-p forces")
print("      partial_x X=-1")

# 9. Hamiltonian preservation gives
# P**2=p**2. The only globally
# smooth solutions compatible with
# diffeomorphism structure are
# P=p and P=-p.
#
# Verify the algebraic factorization.
Pvar = sp.symbols("Pvar",real=True)

assert eq(
    Pvar**2-p**2,
    (Pvar-p)*(Pvar+p)
)

print("PASS: momentum branches")
print("      factor exactly")

# 10. The additive Hamiltonian
# constant must vanish globally.
#
# If H composed T = H+C and T is
# surjective, both functions have
# the same range [m,infinity).
# Translation of this range by C
# is equal only for C=0.
#
# Check the minimum condition
# symbolically.

C = sp.symbols("C",real=True)

assert sp.solve(
    sp.Eq(m+C,m),C
) == [0]

print("PASS: global Hamiltonian")
print("      range excludes nonzero")
print("      additive constants")

# 11. Test concrete nonpolynomial g.
g_example = sp.sin(p)+p**3

for eps in (1,-1):
    T = sp.Matrix([
        eps*x+g_example,
        eps*p
    ])

    assert is_symplectic(T)
    assert commutes(T)
    assert eq(H.subs(p,T[1]),H)

print("PASS: nonlinear examples")
print("      satisfy all conditions")

# 12. Spatial translation is included.
T_translation = sp.Matrix([x+a,p])

assert is_symplectic(T_translation)
assert commutes(T_translation)

print("PASS: spatial translations")
print("      belong to classified family")

# 13. Time translation along the
# free flow is included by g=t0*v.
t0 = sp.symbols("t0",real=True)

T_time = sp.Matrix([
    x+t0*v,p
])

assert is_symplectic(T_time)
assert commutes(T_time)

print("PASS: free-flow time shifts")
print("      belong to classified family")

# 14. Canonical reciprocal scaling
# is generally outside the family.
scale = sp.symbols("scale",positive=True)

T_scale = sp.Matrix([
    scale*x,p/scale
])

assert is_symplectic(T_scale)

assert not eq(
    apply(T_scale,flow(sp.Matrix([x,p])))
        .subs({scale:2,p:1,m:1,t:1}),
    flow(T_scale)
        .subs({scale:2,p:1,m:1,t:1})
)

print("PASS: canonical scaling")
print("      is not generally a symmetry")

# 15. Momentum translation is
# generally outside the family.
T_shift = sp.Matrix([x,p+b])

assert is_symplectic(T_shift)

assert not eq(
    apply(T_shift,flow(sp.Matrix([x,p])))
        .subs({b:1,p:0,m:1,t:1}),
    flow(T_shift)
        .subs({b:1,p:0,m:1,t:1})
)

print("PASS: momentum translations")
print("      are not generally symmetries")

# 16. Composition of two plus-branch
# transformations remains plus-branch.
g1 = sp.Function("g1")(p)
g2 = sp.Function("g2")(p)

A = sp.Matrix([x+g1,p])
B = sp.Matrix([x+g2,p])

composed = apply(A,B)

assert eq(
    composed,
    sp.Matrix([x+g1+g2,p])
)

print("PASS: plus-branch transformations")
print("      are closed under composition")

# 17. The minus branch squared is
# a plus-branch transformation.
R = sp.Matrix([-x+g,-p])

R2 = apply(R,R)

assert eq(
    R2,
    sp.Matrix([
        x-g+g.subs(p,-p),
        p
    ])
)

print("PASS: two reversal-branch")
print("      transformations compose")
print("      into the plus branch")

print("\n"+"="*64)
print("CHECKPOINT 097: ALL TESTS PASSED")
print("="*64)

print("MATHEMATICAL RESULT:")
print("- Under the stated global smoothness")
print("  and symplectic assumptions, the")
print("  free-flow centralizer consists of")
print("  T(x,p)=(eps*x+g(p),eps*p).")
print("- eps is either +1 or -1.")
print("- g is an arbitrary smooth function.")
print("- Both branches preserve H and omega.")
print("- Both branches commute with free flow.")
print("- Spatial and temporal translations")
print("  are included as special cases.")

print("\nPROOF DEPENDENCIES:")
print("- A global symplectomorphism is assumed.")
print("- The complete free flow is fixed.")
print("- Commutation implies H composed T")
print("  differs from H by a constant.")
print("- Surjectivity forces that constant")
print("  to vanish by Hamiltonian range.")
print("- Smoothness forces a consistent")
print("  momentum-sign branch.")
print("- Symplecticity determines X.")

print("\nNOT ESTABLISHED:")
print("- This group is uniquely selected")
print("  by a deeper physical principle.")
print("- Interactions or new field dynamics.")
print("- A noncommutative observable algebra.")
print("- Quantum measurement or Born rule.")
print("- A preferred action scale.")
print("- Emergent quantum mechanics.")
