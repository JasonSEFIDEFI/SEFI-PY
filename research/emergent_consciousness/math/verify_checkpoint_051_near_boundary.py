"""Checkpoint 051: exact near-boundary rational certificate for midpoint.

This certifies a lower bound, not equality with the spectral upper bound.
Author: Jason Duran Dutton.
"""
import itertools
import sympy as sp

A=sp.Matrix([[sp.Rational(3,5),-sp.Rational(2,5)],[0,sp.Rational(1,2)]])
P=sp.Matrix([[sp.Rational(81,341),-sp.Rational(125,874)],
             [-sp.Rational(125,874),sp.Rational(260,341)]])
radius=sp.Rational(76923,500000) # 0.153846
upper=sp.Rational(2,13)
assert radius < upper
assert P[0,0]>0 and sp.factor(P.det())>0
for signs in itertools.product((-1,1),repeat=4):
    M=A+radius*sp.Matrix(2,2,signs)
    H=P-M.T*P*M
    assert H[0,0]>0 and H[1,1]>0 and sp.factor(H.det())>0, signs
print("PASS: P positive definite, exact")
print("PASS: all 16 vertex margins positive definite, exact")
print("CERTIFIED: arbitrary-switching stability for 0<=epsilon<=76923/500000")
print("UPPER: epsilon_critical <= 2/13 from exact spectral witness")
print("UNRESOLVED: endpoint epsilon=2/13 and sharpness")
