import itertools
import numpy as np
import cvxpy as cp
import sympy as sp

print("Checkpoint 038: Boundary metric feasibility")

e = sp.Rational(7, 145)

S = sp.Matrix([[1, 4], [0, 1]])
B = S * sp.diag(sp.Rational(4,5),
                sp.Rational(3,10)) * S.inv()
C = S * sp.diag(sp.Rational(2,5),
                sp.Rational(7,10)) * S.inv()

vertices = []

for name, base in [("B", B), ("C", C)]:
    for signs in itertools.product((-1,1), repeat=4):
        M = base + e * sp.Matrix(2,2,signs)
        vertices.append((name, signs, np.array(M,dtype=float)))

P = cp.Variable((2,2), symmetric=True)
delta = cp.Variable()

constraints = [
    cp.trace(P) == 1,
    P - delta*np.eye(2) >> 0
]

for _,_,M in vertices:
    constraints.append(P-M.T@P@M >> 0)

problem = cp.Problem(cp.Maximize(delta),constraints)

problem.solve(
    solver=cp.CLARABEL,
    tol_gap_abs=1e-11,
    tol_gap_rel=1e-11,
    tol_feas=1e-11,
    max_iter=500
)

print("solver_status =", problem.status)
print("optimized_delta =", delta.value)

if P.value is None:
    print("No boundary candidate returned")
    print("Checkpoint 038 remains unresolved")
    raise SystemExit

Q = (P.value+P.value.T)/2

print("candidate_P =")
print(Q)
print("P_eigenvalues =",np.linalg.eigvalsh(Q))

results = []

for name,signs,M in vertices:
    H = Q-M.T@Q@M
    margin = np.linalg.eigvalsh(H)[0]
    results.append((margin,name,signs))

results.sort(key=lambda x:x[0])

print("Five smallest boundary margins:")
for margin,name,signs in results[:5]:
    print(name,signs,f"{margin:.12e}")

print("Checkpoint 038: numerical diagnosis complete")
print("No exact boundary certificate claimed")
