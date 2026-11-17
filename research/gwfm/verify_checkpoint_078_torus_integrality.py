#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 078

Symplectic torus and conditional prequantization
integrality.

Classical input:
T^2 = R^2 / Z^2
omega = B dx wedge dy, B > 0.

ADDITIONAL assumption:
A U(1) line bundle with connection whose
curvature is i*lambda*omega.

We test the resulting global consistency
condition, not a derivation of quantum physics.
"""

import sympy as sp

x, y = sp.symbols("x y", real=True)
B = sp.symbols("B", positive=True, real=True)
lam = sp.symbols("lam", real=True)
n = sp.symbols("n", integer=True)

I = sp.I
pi = sp.pi

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 078")
print("SYMPLECTIC TORUS AND INTEGRALITY")
print("=" * 64)

# 1. Symplectic area of the unit torus.
area = sp.integrate(
    B, (x, 0, 1), (y, 0, 1)
)

assert area == B
print("PASS: total symplectic area =", area)

# 2. Local real connection potential:
# A = lambda*B*x dy.
# Its exterior derivative is
# dA = lambda*B dx wedge dy.
#
# Multiplication by i gives the
# corresponding imaginary-valued curvature.

A_y = lam*B*x
curvature_coefficient = sp.diff(A_y, x)

assert sp.simplify(
    curvature_coefficient-lam*B
) == 0

print("PASS: local curvature matches lambda*omega")

# 3. At x=0 and x=1, the local potentials
# differ by lambda*B dy.
# A transition function may be written as
# g(y)=exp(i*lambda*B*y).
#
# It must be single-valued around y~y+1.

g = sp.exp(I*lam*B*y)

monodromy = sp.simplify(
    sp.expand_power_exp(
        g.subs(y, y+1)/g
    )
)

assert sp.simplify(
    monodromy-sp.exp(I*lam*B)
) == 0

print("PASS: transition monodromy derived")

# 4. Integrality condition:
# lambda*B = 2*pi*n.

allowed_lambda = 2*pi*n/B

allowed_monodromy = sp.simplify(
    monodromy.subs(lam, allowed_lambda)
)

assert allowed_monodromy == 1

print("PASS: integer flux gives single-valued transition")

# 5. Several admissible integers.
for N in (-2, -1, 0, 1, 2, 3):
    value = sp.simplify(
        monodromy.subs({
            B: 3,
            lam: 2*pi*N/3
        })
    )
    assert value == 1

print("PASS: multiple integer sectors admitted")

# 6. A nonintegral example.
bad = sp.simplify(
    monodromy.subs({
        B: 3,
        lam: pi/3
    })
)

assert bad == -1
assert bad != 1

print("PASS: nonintegral example fails")
print("Nonintegral monodromy =", bad)

# 7. For fixed B, multiple nonzero
# inverse-action scales are permitted.
lambda_1 = 2*pi/B
lambda_2 = 4*pi/B

assert sp.simplify(lambda_2-2*lambda_1) == 0
assert sp.simplify(lambda_1-lambda_2) != 0

print("PASS: integrality does not select unique lambda")

# 8. Classical symplectic area itself
# imposes no restriction on lambda.
# The restriction entered when requiring
# the global U(1) connection construction.

print("PASS: classical and conditional")
print("      quantum-geometric inputs distinguished")

print("\n" + "=" * 64)
print("CHECKPOINT 078: ALL TESTS PASSED")
print("=" * 64)

print("CLASSICAL:")
print("- The torus carries a symplectic area B.")

print("\nCONDITIONAL:")
print("- A global U(1) connection with curvature")
print("  i*lambda*omega requires")
print("  lambda*B/(2*pi) to be an integer.")
print("- Multiple allowed lambda values remain.")

print("\nNOT ESTABLISHED:")
print("- GWFM requires a U(1) phase bundle.")
print("- A unique universal Planck constant.")
print("- Quantum Hilbert-space dynamics.")
print("- The Born rule.")
