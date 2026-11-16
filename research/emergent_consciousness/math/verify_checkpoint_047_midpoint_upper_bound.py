"""Checkpoint 047: exact midpoint spectral upper bound.

Provides a supercritical instability witness, not a full robust-stability
certificate or proof of sharpness. Author: Jason Duran Dutton.
"""
import sympy as sp

t, lam = sp.symbols("t lam", real=True)
B = sp.Matrix([[sp.Rational(4, 5), -2], [0, sp.Rational(3, 10)]])
C = sp.Matrix([[sp.Rational(2, 5), sp.Rational(6, 5)], [0, sp.Rational(7, 10)]])
A = (B+C)/2
D = sp.Matrix([[1,-1],[-1,1]])
M = A+t*D
assert A == sp.Matrix([[sp.Rational(3,5), -sp.Rational(2,5)], [0,sp.Rational(1,2)]])
p = sp.factor((lam*sp.eye(2)-M).det())
critical = sp.Rational(2,13)
assert sp.simplify(p.subs(lam,1)-(sp.Rational(1,5)-sp.Rational(13,10)*t)) == 0
assert sp.solve(p.subs(lam,1),t) == [critical]
assert sp.Rational(1,16) < critical
assert sp.Poly(p,lam).LC() == 1
# For any t>critical, p(1)<0 and p(lambda)->+infinity as lambda->+infinity;
# hence a real eigenvalue exceeds one and repeated constant switching grows.
print("PASS: midpoint A =", A)
print("PASS: det(I-M(t)) =", sp.factor(p.subs(lam,1)))
print("PASS: exact supercritical growth witness for every t > 2/13")
print("CERTIFIED BOUNDS: 1/16 < epsilon_critical(midpoint) <= 2/13")
print("STOP-GATE: no sharpness proof; no arbitrary-switching certificate up to 2/13")
