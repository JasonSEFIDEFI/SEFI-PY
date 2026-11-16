"""Checkpoint 052: no positive-definite common quadratic endpoint metric.

This does NOT establish endpoint instability or unbounded switching.
Author: Jason Duran Dutton.
"""
import sympy as sp

t=sp.Rational(2,13)
A=sp.Matrix([[sp.Rational(3,5),-sp.Rational(2,5)],
             [0,sp.Rational(1,2)]])
D=sp.Matrix([[1,-1],[-1,1]])
M=A+t*D
v=sp.Matrix([8,5])
assert M*v==v
assert sp.factor(M.det())==sp.Rational(53,130)
p,q,r=sp.symbols('p q r',real=True)
P=sp.Matrix([[p,q],[q,r]])
H=P-M.T*P*M
# If H is positive semidefinite and v^T H v = 0, then H v = 0.
assert sp.expand((v.T*H*v)[0])==0
equations=list(H*v)
solution=sp.solve(equations,[p,q],dict=True)
assert len(solution)==1
assert sp.simplify(solution[0][p]-sp.Rational(25,64)*r)==0
assert sp.simplify(solution[0][q]+sp.Rational(5,8)*r)==0
forced=P.subs(solution[0])
assert sp.factor(forced.det())==0
print('PASS: exact endpoint matrix has eigenvector (8,5) at eigenvalue 1')
print('PASS: PSD kernel identity forces H*v=0')
print('PASS: resulting P must equal',forced)
print('PROVED: no positive-definite P can certify even this endpoint vertex')
print('LIMIT: endpoint boundedness and general nonquadratic certificates unresolved')
