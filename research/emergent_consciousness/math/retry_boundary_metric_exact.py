import itertools
import sympy as sp

print("Checkpoint 039: Corrected exact search")

e = sp.Rational(7, 145)
S = sp.Matrix([[1, 4], [0, 1]])

B = S * sp.diag(
    sp.Rational(4, 5),
    sp.Rational(3, 10)
) * S.inv()

C = S * sp.diag(
    sp.Rational(2, 5),
    sp.Rational(7, 10)
) * S.inv()

vertices = []

for name, base in (("B", B), ("C", C)):
    for signs in itertools.product((-1, 1), repeat=4):
        M = base + e * sp.Matrix(2, 2, signs)
        vertices.append((name, signs, M))

def positive_definite(A):
    return A[0, 0] > 0 and A.det() > 0

def positive_semidefinite(A):
    return (
        A[0, 0] >= 0
        and A[1, 1] >= 0
        and A.det() >= 0
    )

found = None

# Search rational r near the numerical candidate.
for k in range(9400, 9701):
    r = sp.Rational(k, 10000)
    q = (608*r - 594) / 145

    P = sp.Matrix([
        [1-r, q],
        [q, r]
    ])

    if not positive_definite(P):
        continue

    if all(
        positive_semidefinite(P-M.T*P*M)
        for _, _, M in vertices
    ):
        found = (r, P)
        break

if found is None:
    print("NO EXACT CERTIFICATE FOUND IN SEARCH")
    print("Checkpoint 039 remains unresolved")
else:
    r, P = found

    print("r =", r)
    print("Exact P =")
    sp.pprint(P)

    assert positive_definite(P)

    zero_determinant_vertices = []

    for name, signs, M in vertices:
        H = P-M.T*P*M

        assert positive_semidefinite(H)

        if H.det() == 0:
            zero_determinant_vertices.append(
                (name, signs)
            )

    print("exact_vertices_verified =", len(vertices))
    print("singular_margin_vertices =",
          zero_determinant_vertices)

    print("PASS: exact positive-definite P")
    print("PASS: all 32 exact semidefinite inequalities")
    print("PASS: Checkpoint 039 verified")
