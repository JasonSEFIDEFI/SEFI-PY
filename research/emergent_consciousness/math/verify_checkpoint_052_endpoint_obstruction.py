"""Checkpoint 052 corrected: endpoint spectral kernel and single-vertex PSD.

Correction: original script used an incorrect eigenvector and incorrectly
claimed singularity of P. This version proves neither existence nor
nonexistence of a COMMON endpoint certificate across all vertices.
Author: Jason Duran Dutton.
"""
import sympy as sp

t=sp.Rational(2,13)
A=sp.Matrix([[sp.Rational(3,5),-sp.Rational(2,5)],
             [0,sp.Rational(1,2)]])
D=sp.Matrix([[1,-1],[-1,1]])
M=A+t*D
v=sp.Matrix([-9,4])
assert M*v==v
assert sp.factor(M.det())==sp.Rational(53,130)
p,q,r=sp.symbols('p q r',real=True)
P=sp.Matrix([[p,q],[q,r]])
H=P-M.T*P*M
assert sp.expand((v.T*H*v)[0])==0
# PSD(H) implies H*v=0. The kernel equation gives ONE linear
# condition on P, not a forced singularity.
assert all(sp.simplify(z.subs(p,-13*q/sp.Integer(72)+5*r/sp.Integer(18)))==0
           for z in H*v)
Pexample=sp.diag(sp.Rational(5,18),1)
Hexample=Pexample-M.T*Pexample*M
assert Pexample.det()>0 and Pexample[0,0]>0
assert Hexample[0,0]>0 and Hexample[1,1]>0
assert sp.factor(Hexample.det())==0
assert Hexample*v==sp.zeros(2,1)
print("PASS: corrected unit eigenvector (-9,4) and other eigenvalue 53/130")
print("PASS: kernel constraint p = -13*q/72 + 5*r/18")
print("PASS: positive-definite P=diag(5/18,1) yields PSD margin for critical vertex")
print("CORRECTION: previous claim of forced singular P was false")
print("UNRESOLVED: common quadratic feasibility across all 16 endpoint vertices")
print("UNRESOLVED: endpoint uniform boundedness under arbitrary switching")
