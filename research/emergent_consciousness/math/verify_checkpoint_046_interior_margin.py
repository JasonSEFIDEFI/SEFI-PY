"""Checkpoint 046: strict interior interpolation and an exact larger-radius witness.

Author: Jason Duran Dutton. All computational tests use exact rational
arithmetic. The conclusion for every fixed interior s also uses
compactness and continuity; this script does not find exact thresholds.
"""
import itertools
import sympy as sp

B = sp.Matrix([[sp.Rational(4, 5), -2], [0, sp.Rational(3, 10)]])
C = sp.Matrix([[sp.Rational(2, 5), sp.Rational(6, 5)], [0, sp.Rational(7, 10)]])
P = sp.Matrix([[sp.Rational(3, 50), -sp.Rational(562, 3625)],
               [-sp.Rational(562, 3625), sp.Rational(47, 50)]])
e0 = sp.Rational(7, 145)
Q = (B-C).T*P*(B-C)
assert Q[0, 0] > 0 and sp.factor(Q.det()) > 0
assert sp.factor((B-C).det()) != 0

def pd2(H):
    return sp.factor(H[0, 0]) > 0 and sp.factor(H.det()) > 0

def psd2(H):
    return (sp.factor(H[0, 0]) >= 0 and sp.factor(H[1, 1]) >= 0
            and sp.factor(H.det()) >= 0)

for base in (B, C):
    for signs in itertools.product((-1, 1), repeat=4):
        M = base + e0*sp.Matrix(2, 2, signs)
        assert psd2(P-M.T*P*M)

# Every fixed interior s has a strict boundary margin: write A(s)+E
# as the interpolation of B+E and C+E; the exact convexity identity
# contributes s(1-s)*Q, positive definite for 0<s<1.
s = sp.symbols('s', real=True)
E = sp.Matrix(sp.symbols('e0:4')).reshape(2, 2)
X, Y = B+E, C+E
Z = (1-s)*X+s*Y
H = P-Z.T*P*Z
R = ((1-s)*(P-X.T*P*X) + s*(P-Y.T*P*Y)
     + s*(1-s)*Q)
assert (H-R).applyfunc(sp.simplify) == sp.zeros(2)

# Concrete radius larger than the global threshold, at s=1/2.
mid = (B+C)/2
higher = sp.Rational(1, 16)
assert higher > e0
for signs in itertools.product((-1, 1), repeat=4):
    M = mid+higher*sp.Matrix(2, 2, signs)
    assert pd2(P-M.T*P*M), signs

print('PASS: (B-C)^T P (B-C) is exactly positive definite')
print('PASS: convexity identity and 32 endpoint vertex PSD tests')
print('PASS: at s=1/2, all 16 vertices strictly contract for epsilon=1/16 > 7/145')
print('CONCLUSION: each fixed 0<s<1 admits some radius >7/145 by compactness and continuity')
print('LIMIT: no exact individual critical radius has been established')
