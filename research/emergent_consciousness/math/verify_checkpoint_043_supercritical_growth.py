"""Checkpoint 043: exact supercritical growing-mode witness.

Author: Jason Duran Dutton
No physical interpretation is inferred.
"""
import sympy as sp

t, lam = sp.symbols("t lam", real=True)
B = sp.Matrix([[sp.Rational(4,5), -2], [0, sp.Rational(3,10)]])
D = sp.Matrix([[1,-1],[-1,1]])
M = B + t*D
char = sp.factor((lam*sp.eye(2)-M).det())
assert sp.factor(char.subs(lam,1)) == (7-145*t)/50
assert sp.Poly(char,lam).LC() == 1
assert sp.solve(char.subs(lam,1),t) == [sp.Rational(7,145)]
print("PASS: det(I-M(t))=(7-145*t)/50")
print("For any t>7/145, characteristic polynomial at 1 is negative")
print("and tends to +infinity as real lambda tends to +infinity.")
print("Thus M(t) has a real eigenvalue lambda>1 and an exponentially")
print("growing trajectory under constant admissible switching.")
print("Characteristic polynomial:", char)
