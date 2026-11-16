"""Checkpoint 054: exact sharp midpoint stability threshold.

Author: Jason Duran Dutton.

Family F(eps) = {A+E: |E_ij| <= eps}, arbitrary switching.
Proves uniform exponential stability iff 0 <= eps < 2/13,
and uniform boundedness at eps=2/13.
Relies on exact 16-vertex endpoint certificate, independently
rechecked here, and strict contraction of nominal A.
"""
import itertools
import sympy as sp

A = sp.Matrix([[sp.Rational(3,5), -sp.Rational(2,5)],
               [0, sp.Rational(1,2)]])
P = sp.Matrix([[sp.Rational(239,720), -sp.Rational(3,10)],
               [-sp.Rational(3,10), 1]])
eps_star = sp.Rational(2,13)
assert P[0,0] > 0 and P.det() == sp.Rational(871,3600)

H0 = sp.simplify(P - A.T * P * A)
assert H0[0,0] == sp.Rational(239,1125)
assert H0.det() == sp.Rational(342047,3240000) > 0
assert H0[0,0] > 0  # H0 positive definite: ||A||_P < 1

for signs in itertools.product((-1,1),repeat=4):
    M = A + eps_star * sp.Matrix(2,2,signs)
    H = sp.simplify(P - M.T*P*M)
    assert H[0,0] >= 0 and H[1,1] >= 0 and H.det() >= 0, signs

# For every 0 <= eps < eps_star and each E with |E_ij|<=eps,
# write theta=eps/eps_star and A+E=(1-theta)*A+theta*(A+E/theta).
# For theta>0, A+E/theta lies in the endpoint box.
# The induced P-operator norm is convex, so
# ||A+E||_P <= (1-theta)*||A||_P + theta*1 < 1.
# The bound is uniform in E and switching sequence.
# At eps=0, A alone strictly contracts.

D = sp.Matrix([[1,-1],[-1,1]])
Mstar = A + eps_star*D
v = sp.Matrix([-9,4])
assert Mstar*v == v
assert sp.factor(Mstar.det()) == sp.Rational(53,130)

# For eps >= eps_star, Mstar is still admissible (nested boxes),
# so its eigenvalue one rules out uniform exponential stability.
# At eps_star, common P bounds arbitrary switching products.

print("PASS: P positive definite, determinant =",P.det())
print("PASS: nominal A strictly P-contracts, H0 determinant =",H0.det())
print("PASS: all 16 endpoint vertices P-nonexpansive (exact arithmetic)")
print("THEOREM: uniform exponential stability iff 0 <= epsilon < 2/13")
print("THEOREM: uniform boundedness at epsilon = 2/13")
print("OBSTRUCTION: admissible unit eigenvalue for all epsilon >= 2/13")
print("SCOPE: specified midpoint family only; no physical SEFI claim")
