"""Checkpoint 042: necessity and completeness of trace-normalized boundary metrics.

Exact linear algebra only. This is a common quadratic certificate result,
not a claim about nonquadratic certificates or physical SEFI phenomena.
Author: Jason Duran Dutton
"""
import sympy as sp

e = sp.Rational(7, 145)
B = sp.Matrix([[sp.Rational(4, 5), -2], [0, sp.Rational(3, 10)]])
D = sp.Matrix([[1, -1], [-1, 1]])
M = B + e*D
v = (M-sp.eye(2)).nullspace()
assert len(v) == 1
v = v[0]
assert M*v == v

p, q, r = sp.symbols("p q r", real=True)
P = sp.Matrix([[p, q], [q, r]])
H = P-M.T*P*M
# If H is PSD and Mv=v, then v.T*H*v=0 forces H*v=0.
assert sp.simplify((v.T*H*v)[0]) == 0
sol = sp.linsolve(list(H*v)+[p+r-1], (p,q,r))
expected = (1-r, (608*r-594)/145, r)
assert sol == sp.FiniteSet(expected), sol
P_family = sp.Matrix([[expected[0],expected[1]],[expected[1],expected[2]]])
assert sp.simplify((P_family-M.T*P_family*M)*v) == sp.zeros(2,1)
assert sp.simplify(P_family.trace()-1) == 0
print("PASS: every positive-definite common quadratic boundary metric")
print("must, after positive trace normalization, lie in P(r).")
print("P(r) =", P_family)
print("Combine with checkpoint 041 exact vertex inequalities for full interval.")
