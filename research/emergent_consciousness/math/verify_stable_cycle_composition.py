import sympy as sp
import numpy as np

print("Checkpoint 021: Stable-cycle composition")

u, v = sp.symbols("u v", positive=True)

def cycle(a, b):
    return sp.exp(-(a+b))*sp.Matrix([
        [1, 10*a],
        [10*b, 1+100*a*b]
    ])

small = sp.Rational(1,100)
large = sp.Integer(1)

B = cycle(small,large)
C = cycle(large,small)
T = C*B

# Exact structural identities.
assert sp.simplify(C-B.T) == sp.zeros(2)
assert sp.simplify(B.det()-sp.exp(-sp.Rational(202,100))) == 0

# Exact eigenvalue formula for a 2x2 matrix.
def eigenvalues_2x2(M):
    tr = sp.trace(M)
    det = M.det()
    disc = sp.simplify(tr**2-4*det)
    return (
        sp.N((tr+sp.sqrt(disc))/2, 40),
        sp.N((tr-sp.sqrt(disc))/2, 40)
    )

eig_B = eigenvalues_2x2(B)
eig_C = eigenvalues_2x2(C)
eig_T = eigenvalues_2x2(T)

assert 0 < eig_B[1] < eig_B[0] < 1
assert 0 < eig_C[1] < eig_C[0] < 1
assert eig_T[0] > 1

print("PASS: individual cycles asymptotically stable")
print("B eigenvalues:", [float(x) for x in eig_B])
print("C eigenvalues:", [float(x) for x in eig_C])
print("T eigenvalues:", [float(x) for x in eig_T])

# Independent numerical matrix calculation.
Bn = np.array(B.evalf(17), dtype=float)
Cn = np.array(C.evalf(17), dtype=float)
Tn = Cn @ Bn

numerical_eigs = np.linalg.eigvals(Tn)

assert np.isclose(
    max(abs(numerical_eigs)),
    float(eig_T[0]),
    rtol=1e-12
)

print("PASS: independent matrix-product verification")

# Test the expanding eigenvector of the product.
values, vectors = np.linalg.eig(Tn)
index = np.argmax(abs(values))
x = np.real(vectors[:,index])
x = x/np.linalg.norm(x)

initial_norm = np.linalg.norm(x)
growth = float(eig_T[0])

for n in range(1,9):
    x = Tn @ x
    actual = np.linalg.norm(x)
    expected = growth**n*initial_norm

    assert np.isclose(
        actual,expected,rtol=1e-10
    )

    print(
        f"pair={n:2d} "
        f"norm={actual:.6e} "
        f"expected={expected:.6e}"
    )

print("PASS: alternating stable cycles are unstable")
print("PASS: Checkpoint 021 verified")
