"""Checkpoint 053: exact endpoint common quadratic feasibility.

Author: Jason Duran Dutton.
No floating-point arithmetic or solver conclusions enter this certificate.
"""
import itertools
import sympy as sp

A = sp.Matrix([[sp.Rational(3, 5), -sp.Rational(2, 5)],
               [sp.Integer(0), sp.Rational(1, 2)]])
epsilon = sp.Rational(2, 13)
P = sp.Matrix([[sp.Rational(239, 720), -sp.Rational(3, 10)],
               [-sp.Rational(3, 10), sp.Integer(1)]])
assert P[0, 0] > 0 and P.det() > 0, "P must be positive definite"

critical = (1, -1, -1, 1)
zero_margins = []
for signs in itertools.product((-1, 1), repeat=4):
    D = sp.Matrix(2, 2, signs)
    M = A + epsilon * D
    H = sp.simplify(P - M.T * P * M)
    assert H == H.T
    assert H[0, 0] >= 0 and H[1, 1] >= 0 and H.det() >= 0, (
        f"PSD margin failed for {signs}: {H}"
    )
    if H.det() == 0:
        zero_margins.append(signs)

assert zero_margins == [critical], zero_margins
Mcritical = A + epsilon * sp.Matrix(2, 2, critical)
v = sp.Matrix([-9, 4])
assert Mcritical * v == v
assert (P - Mcritical.T * P * Mcritical) * v == sp.zeros(2, 1)

# Any entrywise uncertainty matrix is a convex combination of its 16 vertices.
# Since ||.||_P is a norm on matrices, convexity gives ||M||_P <= 1
# throughout the box; arbitrary products therefore satisfy ||product||_P <= 1.
print("PASS: P positive definite, exact determinant =", P.det())
print("PASS: all 16 endpoint vertex Lyapunov margins positive semidefinite")
print("PASS: exactly one vertex margin is singular")
print("PASS: endpoint family is uniformly bounded under arbitrary switching")
print("LIMIT: endpoint is not uniformly exponentially stable (unit eigenvalue)")
print("LIMIT: this result alone does not prove sharp subcritical stability up to 2/13")
print("LIMIT: no physical SEFI correspondence follows from this certificate")
