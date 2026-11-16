"""Checkpoint 045: exact common quadratic certificate for convex interpolation.

Proves a sufficient boundary certificate for A(s)+E, 0<=s<=1,
with |E_ij|<=7/145. Does not assert a sharp threshold for each s.
Author: Jason Duran Dutton
"""
import itertools
import sympy as sp

eps = sp.Rational(7, 145)
B = sp.Matrix([[sp.Rational(4, 5), -2], [0, sp.Rational(3, 10)]])
C = sp.Matrix([[sp.Rational(2, 5), sp.Rational(6, 5)], [0, sp.Rational(7, 10)]])
P = sp.Matrix([[sp.Rational(3,50), -sp.Rational(562,3625)],
               [-sp.Rational(562,3625), sp.Rational(47,50)]])
assert P[0,0] > 0 and P.det() > 0

def psd2(H):
    H = H.applyfunc(sp.factor)
    assert H == H.T
    return H[0,0] >= 0 and H[1,1] >= 0 and sp.factor(H.det()) >= 0

vertices = []
for name, A in (("B",B), ("C",C)):
    for signs in itertools.product((-1,1), repeat=4):
        M = A + eps*sp.Matrix(2,2,signs)
        H = P - M.T*P*M
        assert psd2(H), (name,signs)
        vertices.append(M)
assert len(vertices) == 32

s = sp.symbols("s", real=True)
X = sp.Matrix(sp.symbols("x0:4")).reshape(2,2)
Y = sp.Matrix(sp.symbols("y0:4")).reshape(2,2)
Z = (1-s)*X+s*Y
lhs = P-Z.T*P*Z
rhs = ((1-s)*(P-X.T*P*X) + s*(P-Y.T*P*Y)
       +s*(1-s)*(X-Y).T*P*(X-Y))
assert (lhs-rhs).applyfunc(sp.simplify) == sp.zeros(2)

# Any A(s)+E is a convex combination of the 32 verified vertices:
# E itself is a convex combination of the 16 box vertices, and
# A(s)+E=(1-s)(B+E)+s(C+E).
# Matrix convexity then establishes the certificate for the hull.
print("PASS: P positive definite (exact rational arithmetic)")
print("PASS: 32 endpoint vertex PSD inequalities")
print("PASS: symbolic matrix convexity identity")
print("CONCLUSION: P-M.T*P*M PSD for all A(s)+E, 0<=s<=1,")
print("            |E_ij|<=7/145, including time-varying switching")
print("LIMIT: no sharp interior-s threshold or novel theorem claimed")
