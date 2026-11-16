"""Checkpoint 044: exact simple unit-eigenvalue crossing.

This verifies symbolic identities, not robust stability for a full uncertainty box.
Author: Jason Duran Dutton
"""
import sympy as sp

a, b, d, t, lam = sp.symbols('a b d t lam', real=True)
A = sp.Matrix([[a, b], [0, d]])
D = sp.Matrix([[1, -1], [-1, 1]])
M = A + t * D
delta = 2 - a - d - b
t_star = (1 - a) * (1 - d) / delta
lam2 = a + d - 1 + 2 * t_star
f = sp.factor((lam * sp.eye(2) - M).det())
assert sp.simplify(f.subs(lam, 1) - ((1-a)*(1-d) - t*delta)) == 0
assert sp.simplify(f.subs({lam: 1, t: t_star})) == 0
assert sp.simplify(M.subs(t, t_star).trace() - (1 + lam2)) == 0
assert sp.simplify(f.subs(t, t_star) - (lam-1)*(lam-lam2)) == 0
slope = sp.factor(-sp.diff(f,t).subs({t:t_star,lam:1}) /
                  sp.diff(f,lam).subs({t:t_star,lam:1}))
assert sp.simplify(slope - delta/(1-lam2)) == 0
B = {a:sp.Rational(4,5), b:-2, d:sp.Rational(3,10)}
assert sp.simplify(t_star.subs(B) - sp.Rational(7,145)) == 0
assert sp.simplify(lam2.subs(B) - sp.Rational(57,290)) == 0
assert slope.subs(B) > 0
print('PASS: exact determinant and eigenvalue identities')
print('PASS: simple crossing slope = delta/(1-lambda_2)')
print('PASS: B reproduces t*=7/145 and lambda_2=57/290')
print('B crossing slope =', sp.factor(slope.subs(B)))
print('Scope: conditional local spectral crossing, not box robustness')
