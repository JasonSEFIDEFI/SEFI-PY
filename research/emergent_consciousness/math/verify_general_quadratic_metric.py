import itertools
import numpy as np
import cvxpy as cp

print("Checkpoint 033: General quadratic metric")

S = np.array([[1.,4.],[0.,1.]])
Sinv = np.linalg.inv(S)

B = S @ np.diag([0.8,0.3]) @ Sinv
C = S @ np.diag([0.4,0.7]) @ Sinv

eps = 0.035

patterns = list(
    itertools.product([-1.,1.],repeat=4)
)

vertices = []

for name,base in [("B",B),("C",C)]:
    for pattern in patterns:
        E = eps*np.array(pattern).reshape(2,2)
        vertices.append((name,pattern,base+E))

P = cp.Variable((2,2),symmetric=True)
t = cp.Variable()

constraints = [
    cp.trace(P) == 1,
    P - t*np.eye(2) >> 0
]

for _,_,M in vertices:
    constraints.append(
        P - M.T@P@M - t*np.eye(2) >> 0
    )

problem = cp.Problem(
    cp.Maximize(t),
    constraints
)

problem.solve(
    solver=cp.CLARABEL,
    verbose=False
)

print("solver_status =", problem.status)
print("optimal_margin =", t.value)

if P.value is None or t.value is None:
    print("NO NUMERICAL CANDIDATE RETURNED")
    print("Checkpoint 033 remains unresolved")
    raise SystemExit(0)

Pnum = np.asarray(P.value)
Pnum = (Pnum+Pnum.T)/2

print("candidate_P =")
print(Pnum)

p_eigenvalues = np.linalg.eigvalsh(Pnum)
print("P_eigenvalues =", p_eigenvalues)

worst_margin = float("inf")
worst_vertex = None

for name,pattern,M in vertices:
    H = Pnum-M.T@Pnum@M
    margin = np.linalg.eigvalsh(H)[0]

    if margin < worst_margin:
        worst_margin = margin
        worst_vertex = (name,pattern)

print("verified_min_P_eigenvalue =", p_eigenvalues[0])
print("verified_worst_margin =", worst_margin)
print("limiting_vertex =", worst_vertex)

tolerance = 1e-8

if (
    p_eigenvalues[0] > tolerance
    and worst_margin > tolerance
):
    print("PASS: strict common quadratic certificate")
    print("PASS: Checkpoint 033 verified")
else:
    print("NO STRICT CERTIFICATE VERIFIED")
    print("Checkpoint 033 remains unresolved")
