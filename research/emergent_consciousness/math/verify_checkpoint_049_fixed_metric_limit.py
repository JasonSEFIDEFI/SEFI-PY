"""Checkpoint 049: exact maximum nonexpansive radius for fixed midpoint P."""
import itertools
import sympy as sp

x = sp.symbols('x', real=True)
B = sp.Matrix([[sp.Rational(4,5), -2], [0, sp.Rational(3,10)]])
C = sp.Matrix([[sp.Rational(2,5), sp.Rational(6,5)], [0, sp.Rational(7,10)]])
A = (B+C)/2
P = sp.Matrix([[sp.Rational(3,50), -sp.Rational(562,3625)],
               [-sp.Rational(562,3625), sp.Rational(47,50)]])
active = (1,1,-1,-1)
critical = (-sp.Integer(10413401)+7*sp.sqrt(597925520449617))/1795734335
assert sp.Rational(89520,1000000) < critical < sp.Rational(89521,1000000)
assert P[0,0]>0 and P.det()>0
checked=0
for signs in itertools.product((-1,1),repeat=4):
    M = A+x*sp.Matrix(2,2,signs)
    H = P-M.T*P*M
    for name, expr in [('h11',H[0,0]), ('h22',H[1,1]), ('det',H.det())]:
        poly = sp.Poly(sp.factor(expr),x)
        assert poly.eval(0)>0, (signs,name,'at zero')
        positive_roots = [z for z in sp.real_roots(poly) if z>0]
        if signs==active and name=='det':
            assert positive_roots==[critical]
        else:
            assert all(z>critical for z in positive_roots), (signs,name,positive_roots)
        checked+=1
# All minors strictly positive before the first positive root; at
# critical, only the active determinant vanishes, so all margins PSD.
active_M=A+x*sp.Matrix(2,2,active)
active_det=sp.factor((P-active_M.T*P*active_M).det())
assert sp.simplify(active_det + (8978671675*x*x+104134010*x-81275696)/5256250000)==0
assert sp.simplify(active_det.subs(x,critical))==0
assert sp.diff(active_det,x).subs(x,critical)<0
print('PASS:',checked,'exact principal-minor polynomials checked')
print('PASS: all 16 vertices strictly contract for 0 <= epsilon < epsilon_P')
print('PASS: all 16 vertices nonexpansive at epsilon_P')
print('PASS: active vertex loses PSD immediately above epsilon_P')
print('epsilon_P =',critical)
print('decimal reference =',critical.evalf(18))
print('LIMIT: exact maximum for THIS fixed P, not the sharp system threshold')
