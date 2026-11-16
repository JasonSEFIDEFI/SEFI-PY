"""Checkpoint 048: exact rational-radius bracket for the existing metric.

The tested metric is fixed; failure does not exclude another metric.
Author: Jason Duran Dutton
"""
import itertools
import sympy as sp

B = sp.Matrix([[sp.Rational(4,5),-2],[0,sp.Rational(3,10)]])
C = sp.Matrix([[sp.Rational(2,5),sp.Rational(6,5)],[0,sp.Rational(7,10)]])
A = (B+C)/2
P = sp.Matrix([[sp.Rational(3,50),-sp.Rational(562,3625)],
               [-sp.Rational(562,3625),sp.Rational(47,50)]])
assert P[0,0]>0 and P.det()>0

def margins(radius):
    for signs in itertools.product((-1,1),repeat=4):
        M=A+radius*sp.Matrix(2,2,signs)
        H=P-M.T*P*M
        yield signs,sp.factor(H[0,0]),sp.factor(H.det())

good=sp.Rational(1,12)
bad=sp.Rational(1,10)
assert all(h11>0 and det>0 for _,h11,det in margins(good))
failures=[(signs,h11,det) for signs,h11,det in margins(bad)
          if not(h11>0 and det>0)]
assert failures
assert sp.Rational(1,16)<good<bad<sp.Rational(2,13)
print("PASS: existing P strictly contracts all 16 vertices at epsilon=1/12")
print("PASS: existing P fails at epsilon=1/10; failing vertices:",len(failures))
print("One failing vertex and exact minors:",failures[0])
print("CERTIFIED: midpoint switching stability for 0<=epsilon<=1/12")
print("LIMIT: failure of this P at 1/10 does not prove instability")
print("LIMIT: sharp midpoint critical radius remains unknown")
