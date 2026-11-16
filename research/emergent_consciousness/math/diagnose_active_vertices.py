import itertools
import numpy as np
import cvxpy as cp

print("Checkpoint 035: Active-vertex diagnosis")

S = np.array([[1., 4.], [0., 1.]])
Si = np.linalg.inv(S)

bases = {
    "B": S @ np.diag([0.8, 0.3]) @ Si,
    "C": S @ np.diag([0.4, 0.7]) @ Si
}

eps = 0.04827
vertices = []

for name, base in bases.items():
    for signs in itertools.product([-1., 1.], repeat=4):
        M = base + eps*np.array(signs).reshape(2, 2)
        vertices.append((name, signs, M))

P = cp.Variable((2, 2), symmetric=True)
t = cp.Variable()

constraints = [
    cp.trace(P) == 1,
    P - t*np.eye(2) >> 0
]

for _, _, M in vertices:
    constraints.append(
        P - M.T @ P @ M - t*np.eye(2) >> 0
    )

problem = cp.Problem(cp.Maximize(t), constraints)
problem.solve(
    solver=cp.CLARABEL,
    tol_gap_abs=1e-11,
    tol_feas=1e-11,
    tol_gap_rel=1e-11,
    max_iter=500
)

print("status =", problem.status)
print("optimized_t =", t.value)

if P.value is None:
    raise RuntimeError("No candidate metric returned")

Q = (P.value + P.value.T)/2

print("P =")
print(Q)
print("P eigenvalues =", np.linalg.eigvalsh(Q))

results = []

for name, signs, M in vertices:
    H = Q - M.T @ Q @ M
    margin = np.linalg.eigvalsh(H)[0]
    results.append((margin, name, signs))

results.sort(key=lambda item: item[0])

print("\nTen smallest vertex margins:")
for margin, name, signs in results[:10]:
    print(
        f"{name} {signs} "
        f"minimum_eigenvalue={margin:.12e}"
    )

print("\nminimum_margin =", results[0][0])
print("trace_P =", np.trace(Q))

assert np.linalg.eigvalsh(Q)[0] > 0
assert results[0][0] > 1e-9

print("PASS: strict numerical certificate at epsilon=0.04827")
print("Checkpoint 035: active-vertex diagnosis complete")
