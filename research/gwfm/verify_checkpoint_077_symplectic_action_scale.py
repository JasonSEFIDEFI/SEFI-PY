#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 077

Does classical symplectic geometry determine
a unique inverse-action scale?

Classical input:
omega = dx wedge dp

Loop action:
A(C) = integral_C p dx

Conditional phase assignment:
F_lambda(C) = exp(i*lambda*A(C))

No quantum phase law is assumed physically.
"""

import sympy as sp

I = sp.I
lam, s = sp.symbols("lambda s", real=True)
r = sp.symbols("r", positive=True, real=True)

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 077")
print("SYMPLECTIC GEOMETRY AND ACTION-SCALE FREEDOM")
print("=" * 64)

def loop_action(vertices):
    """Exact integral of p dx along polygon edges."""
    assert vertices[0] == vertices[-1]
    total = 0
    for (x1,p1),(x2,p2) in zip(
        vertices[:-1], vertices[1:]
    ):
        # Linear interpolation on each edge.
        total += (p1+p2)*(x2-x1)/2
    return sp.simplify(total)

def transform(vertices, matrix):
    return [
        tuple(matrix*sp.Matrix([x,p]))
        for x,p in vertices
    ]

def reverse(vertices):
    return list(reversed(vertices))

def concatenate(A,B):
    assert A[0] == A[-1] == B[0] == B[-1]
    return A[:-1] + B

def phase(action, scale):
    return sp.exp(I*scale*action)

# 1. Closed rectangular loop.
C = [
    (sp.Integer(0),sp.Integer(0)),
    (sp.Integer(2),sp.Integer(0)),
    (sp.Integer(2),sp.Integer(3)),
    (sp.Integer(0),sp.Integer(3)),
    (sp.Integer(0),sp.Integer(0)),
]

A = loop_action(C)
assert A == -6

print("PASS: exact closed-loop action =", A)

# 2. Orientation reversal.
assert sp.simplify(
    loop_action(reverse(C)) + A
) == 0

print("PASS: orientation reversal changes action sign")

# 3. Symplectic shear:
# x' = x+s*p, p'=p
shear = sp.Matrix([[1,s],[0,1]])

Omega = sp.Matrix([[0,1],[-1,0]])

assert sp.simplify(
    shear.T*Omega*shear-Omega
) == sp.zeros(2)

assert sp.simplify(
    loop_action(transform(C,shear))-A
) == 0

print("PASS: symplectic shear preserves loop action")

# 4. Symplectic scaling:
# x'=r*x, p'=p/r
scaling = sp.diag(r,1/r)

assert sp.simplify(
    scaling.T*Omega*scaling-Omega
) == sp.zeros(2)

assert sp.simplify(
    loop_action(transform(C,scaling))-A
) == 0

print("PASS: canonical scaling preserves loop action")

# 5. Loop concatenation.
D = [
    (sp.Integer(0),sp.Integer(0)),
    (sp.Integer(1),sp.Integer(0)),
    (sp.Integer(1),sp.Integer(2)),
    (sp.Integer(0),sp.Integer(2)),
    (sp.Integer(0),sp.Integer(0)),
]

AD = loop_action(D)
combined = concatenate(C,D)

assert AD == -2
assert loop_action(combined) == A+AD

print("PASS: loop action is additive")

# 6. Conditional phase characters.
for k in (sp.Integer(1),sp.Integer(2),
          sp.Rational(1,2)):
    lhs = phase(A+AD,k)
    rhs = phase(A,k)*phase(AD,k)
    assert sp.simplify(
        sp.expand_power_exp(lhs/rhs)
    ) == 1

print("PASS: multiple phase scales preserve composition")

# 7. Distinct scales give distinct phases
# for a fixed action.
F1 = phase(A,sp.pi/12)
F2 = phase(A,sp.pi/6)

assert sp.simplify(F1) == -I
assert sp.simplify(F2) == -1
assert F1 != F2

print("PASS: same geometric action, distinct phases")
print("scale pi/12 ->", F1)
print("scale pi/6  ->", F2)

# 8. Symbolic invariance holds for arbitrary
# lambda, not only a selected scale.
assert sp.simplify(
    phase(loop_action(transform(C,shear)),lam)
    - phase(A,lam)
) == 0

assert sp.simplify(
    phase(loop_action(transform(C,scaling)),lam)
    - phase(A,lam)
) == 0

print("PASS: symplectic invariance leaves lambda free")

# 9. Units:
# integral p dx has dimensions of action.
# lambda must have dimensions inverse action
# for lambda*A to be dimensionless.
# Dimensional consistency does not set a value.

print("PASS: action-scale freedom retained")

print("\n" + "=" * 64)
print("CHECKPOINT 077: ALL TESTS PASSED")
print("=" * 64)

print("ESTABLISHED:")
print("- Closed-loop action is additive.")
print("- Orientation reversal reverses its sign.")
print("- Tested canonical transformations preserve it.")
print("- Symplectic invariance does not fix lambda.")

print("\nNOT ESTABLISHED:")
print("- A physically required complex phase.")
print("- lambda = 1/hbar.")
print("- A universal quantum action scale.")
print("- The Born rule.")
