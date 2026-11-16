import itertools
import sympy as sp

print("Checkpoint 030: Exact robustness boundary")

e = sp.symbols("e", real=True)

S = sp.Matrix([[1,4],[0,1]])
P = sp.Matrix([[1,-4],[-4,41]])

B = S*sp.diag(
    sp.Rational(4,5),
    sp.Rational(3,10)
)*S.inv()

C = S*sp.diag(
    sp.Rational(2,5),
    sp.Rational(7,10)
)*S.inv()

critical = sp.Rational(13,375)

E = e*sp.Matrix([[1,1],[-1,-1]])
H = P-(B+E).T*P*(B+E)

expected = -(375*e-13)*(625*e+63)/100

assert sp.simplify(H.det()-expected) == 0
assert sp.simplify(H.det().subs(e,critical)) == 0

print("PASS: exact determinant factorization")
print("critical =", critical)
print("decimal =", sp.N(critical,15))

# Verify strict positive definiteness for
# all vertices just below the boundary.
below = critical*sp.Rational(99,100)
above = critical*sp.Rational(101,100)

signs = list(itertools.product([-1,1],repeat=4))

def check_vertices(epsilon):
    minimum_det = None
    failures = 0

    for base in (B,C):
        for pattern in signs:
            perturbation = epsilon*sp.Matrix(
                2,2,pattern
            )
            M = base+perturbation
            G = P-M.T*P*M

            d = sp.factor(G.det())

            if G[0,0] <= 0 or d <= 0:
                failures += 1

            if minimum_det is None or d < minimum_det:
                minimum_det = d

    return failures,minimum_det

fail_below,det_below = check_vertices(below)
fail_above,det_above = check_vertices(above)

assert fail_below == 0
assert fail_above > 0

print("PASS: all 32 vertices strict below boundary")
print("failed_vertices_above =", fail_above)

# Verify that the limiting vertex loses
# positive definiteness at the boundary.
Hcrit = H.subs(e,critical)

assert Hcrit.det() == 0
assert Hcrit[0,0] > 0
assert sp.trace(Hcrit) > 0

print("PASS: exact semidefinite boundary")
print("PASS: Checkpoint 030 verified")
