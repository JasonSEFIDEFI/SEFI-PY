"""Checkpoint 041: exact trace-normalized boundary-metric interval.

Uses rational and algebraic SymPy arithmetic; no floating-point sign tests.
Author: Jason Duran Dutton
"""
import itertools
import sympy as sp

r = sp.symbols("r", real=True)
e = sp.Rational(7, 145)
B = sp.Matrix([[sp.Rational(4, 5), -2], [0, sp.Rational(3, 10)]])
C = sp.Matrix([[sp.Rational(2, 5), sp.Rational(6, 5)], [0, sp.Rational(7, 10)]])
q = (608*r - 594)/145
P = sp.Matrix([[1-r, q], [q, r]])

lower_poly = 827071751*r**2 - 1585543261*r + 758855474
upper_poly = 2571639607*r**2 - 4871334677*r + 2303901468
lower = min(sp.solve(lower_poly, r))
upper = max(sp.solve(upper_poly, r))
middle = sp.Rational(95, 100)

assert sp.Rational(923, 1000) < lower < sp.Rational(924, 1000)
assert sp.Rational(981, 1000) < upper < sp.Rational(982, 1000)
assert lower < middle < upper
assert sp.simplify(P.det() + (533*r-484)*(733*r-729)/21025) == 0
assert sp.Rational(484, 533) < lower < upper < sp.Rational(729, 733)

# Sufficiency: all three 2x2 principal-minor conditions for all 32
# uncertainty-box vertices. Exact algebraic root comparisons exclude
# sign changes in the interior; signs at endpoints and middle are exact.
checks = 0
for name, base in (("B", B), ("C", C)):
    for signs in itertools.product((-1, 1), repeat=4):
        M = base + e*sp.Matrix(2, 2, signs)
        H = sp.simplify(P - M.T*P*M)
        for label, expr in (("h11", H[0,0]), ("h22", H[1,1]),
                            ("det", H.det())):
            poly = sp.Poly(sp.factor(expr), r)
            if poly.is_zero:
                checks += 1
                continue
            for root in sp.real_roots(poly):
                assert not (lower < root < upper), (name, signs, label, root)
            for x in (lower, middle, upper):
                assert sp.sign(poly.eval(x)) >= 0, (name, signs, label, x)
            checks += 1

# Necessity: two active B vertices exclude every r outside [lower,upper]
# within the positive-definite domain (484/533,729/733).
for signs, target, endpoint in (
    ((1,1,-1,-1), lower_poly, lower),
    ((1,-1,1,-1), upper_poly, upper),
):
    M = B + e*sp.Matrix(2, 2, signs)
    expr = sp.factor((P-M.T*P*M).det())
    assert sp.rem(sp.Poly(expr,r), sp.Poly(target,r)).is_zero
    other = [x for x in sp.solve(target,r) if x != endpoint][0]
    if endpoint == lower:
        assert upper < other
    else:
        assert other < lower
    assert sp.sign(expr.subs(r,middle)) > 0
    print("active B vertex", signs, "determinant =", expr)

print("PASS:", checks, "exact scalar inequalities across 32 vertices")
print("trace-one admissible r interval:", lower, upper)
print("decimal reference only:", lower.evalf(18), upper.evalf(18))
