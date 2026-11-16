import itertools
import numpy as np
import cvxpy as cp
import sympy as sp

print("Checkpoint 037: Exact rational certificate")

eps = sp.Rational(1931, 40000)

S = sp.Matrix([[1, 4], [0, 1]])
B = S * sp.diag(sp.Rational(4, 5),
                sp.Rational(3, 10)) * S.inv()
C = S * sp.diag(sp.Rational(2, 5),
                sp.Rational(7, 10)) * S.inv()

vertices = []
for name, base in [("B", B), ("C", C)]:
    for signs in itertools.product((-1, 1), repeat=4):
        M = base + eps * sp.Matrix(2, 2, signs)
        vertices.append((name, signs, M))

P = cp.Variable((2, 2), symmetric=True)
t = cp.Variable()

constraints = [
    cp.trace(P) == 1,
    P - t * np.eye(2) >> 0
]

for _, _, M in vertices:
    A = np.array(M, dtype=float)
    constraints.append(
        P - A.T @ P @ A - t * np.eye(2) >> 0
    )

problem = cp.Problem(cp.Maximize(t), constraints)
problem.solve(
    solver=cp.CLARABEL,
    tol_gap_abs=1e-11,
    tol_gap_rel=1e-11,
    tol_feas=1e-11,
    max_iter=500
)

print("solver_status =", problem.status)

if P.value is None:
    raise RuntimeError("No candidate metric returned")

Q = (P.value + P.value.T) / 2

# Rationalize to 12 decimal places.
scale = 10**12
R = sp.Matrix([
    [sp.Rational(round(float(Q[i, j]) * scale), scale)
     for j in range(2)]
    for i in range(2)
])

assert R == R.T

def positive_definite_2x2(A):
    return A[0, 0] > 0 and A.det() > 0

assert positive_definite_2x2(R)

worst = None
for name, signs, M in vertices:
    H = R - M.T * R * M
    assert H == H.T

    if not positive_definite_2x2(H):
        raise AssertionError(
            f"Exact certificate failed: {name}, {signs}"
        )

    determinant = H.det()
    if worst is None or determinant < worst[0]:
        worst = (determinant, name, signs)

print("rational_P =")
sp.pprint(R)

print("exact_vertices_verified =", len(vertices))
print("minimum_determinant =", sp.N(worst[0], 16))
print("minimum_determinant_vertex =", worst[1:])
print("certified_epsilon =", eps)
print("remaining_gap =", sp.Rational(7,145) - eps)

print("PASS: exact positive-definite metric")
print("PASS: all 32 exact Lyapunov inequalities")
print("PASS: Checkpoint 037 verified")
