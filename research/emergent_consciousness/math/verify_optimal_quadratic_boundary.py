import itertools
import numpy as np
import cvxpy as cp

print("Checkpoint 034: Optimal quadratic boundary")

S = np.array([[1.,4.],[0.,1.]])
Si = np.linalg.inv(S)

bases = {
    "B": S @ np.diag([0.8,0.3]) @ Si,
    "C": S @ np.diag([0.4,0.7]) @ Si
}

patterns = list(itertools.product([-1.,1.],repeat=4))

def solve_at(eps):
    vertices = [
        (name,pattern,base+eps*np.array(pattern).reshape(2,2))
        for name,base in bases.items()
        for pattern in patterns
    ]

    P = cp.Variable((2,2),symmetric=True)
    t = cp.Variable()

    constraints = [
        cp.trace(P)==1,
        P-t*np.eye(2) >> 0
    ]

    for _,_,M in vertices:
        constraints.append(
            P-M.T@P@M-t*np.eye(2) >> 0
        )

    problem = cp.Problem(cp.Maximize(t),constraints)
    problem.solve(solver=cp.CLARABEL)

    if P.value is None or t.value is None:
        return problem.status,None,None,None

    Q = np.asarray(P.value)
    Q = (Q+Q.T)/2

    min_p = np.linalg.eigvalsh(Q)[0]
    margins = [
        np.linalg.eigvalsh(Q-M.T@Q@M)[0]
        for _,_,M in vertices
    ]

    return problem.status,float(t.value),float(min_p),float(min(margins))

# Known feasible lower endpoint.
low = 0.035
status,t,min_p,margin = solve_at(low)

assert status == "optimal"
assert min_p > 1e-8 and margin > 1e-8

# Absolute spectral upper endpoint.
high = 0.056

for iteration in range(18):
    mid = (low+high)/2
    status,t,min_p,margin = solve_at(mid)

    feasible = (
        status in ("optimal","optimal_inaccurate")
        and min_p is not None
        and min_p > 1e-7
        and margin > 1e-7
    )

    print(
        f"{iteration+1:02d} "
        f"epsilon={mid:.10f} "
        f"status={status} "
        f"margin={margin} "
        f"certified={feasible}"
    )

    if feasible:
        low = mid
    else:
        high = mid

print()
print(f"verified_feasible_lower={low:.12f}")
print(f"unverified_upper={high:.12f}")
print(f"numerical_bracket_width={high-low:.12e}")
print("Checkpoint 034: numerical search completed")
