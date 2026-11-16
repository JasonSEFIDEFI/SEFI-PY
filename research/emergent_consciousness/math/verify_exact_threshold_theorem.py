import itertools
import sympy as sp

print("Checkpoint 040: Exact robust stability theorem")

ec = sp.Rational(7, 145)
S = sp.Matrix([[1, 4], [0, 1]])

B = S*sp.diag(sp.Rational(4,5),
              sp.Rational(3,10))*S.inv()
C = S*sp.diag(sp.Rational(2,5),
              sp.Rational(7,10))*S.inv()

P = sp.Matrix([
    [sp.Rational(3,50), -sp.Rational(562,3625)],
    [-sp.Rational(562,3625), sp.Rational(47,50)]
])

def pd(A):
    return A[0,0] > 0 and A.det() > 0

def psd(A):
    return (A[0,0] >= 0 and
            A[1,1] >= 0 and
            A.det() >= 0)

assert pd(P)
print("PASS: P positive definite")

for name,A in (("B",B),("C",C)):
    H = P-A.T*P*A
    assert pd(H)
    print(name, "unperturbed margin determinant =", H.det())

print("PASS: both unperturbed operators strictly contract")

singular = []
for name,A in (("B",B),("C",C)):
    for signs in itertools.product((-1,1),repeat=4):
        M = A+ec*sp.Matrix(2,2,signs)
        H = P-M.T*P*M
        assert psd(H)
        if H.det() == 0:
            singular.append((name,signs))

print("PASS: all 32 boundary inequalities exact")
print("singular_vertices =", singular)

assert singular == [("B",(1,-1,-1,1))]

D = sp.Matrix([[1,-1],[-1,1]])
Mstar = B+ec*D
assert (sp.eye(2)-Mstar).det() == 0
assert sp.Integer(1) in Mstar.eigenvals()

print("PASS: boundary spectral obstruction")
print("exact_threshold =", ec)
print("decimal_threshold =", sp.N(ec,17))

print("THEOREM VERIFIED:")
print("For every 0 <= epsilon < 7/145,")
print("the entire uncertainty family is uniformly")
print("exponentially stable under arbitrary switching.")
print("At epsilon = 7/145, strict asymptotic")
print("stability of the entire family fails.")
print("PASS: Checkpoint 040 verified")
