"""Checkpoint 050: numerical search + exact rational certificate.

Numerical SDP output is exploratory. Only exact SymPy rational checks
constitute a verified certificate. Author: Jason Duran Dutton.
"""
import itertools
import sys
import sympy as sp

try:
    import numpy as np
    import cvxpy as cp
except ImportError as exc:
    print("MISSING DEPENDENCY:", exc)
    print("Install in your active venv: python -m pip install numpy cvxpy sympy")
    sys.exit(2)

B = sp.Matrix([[sp.Rational(4,5),-2],[0,sp.Rational(3,10)]])
C = sp.Matrix([[sp.Rational(2,5),sp.Rational(6,5)],[0,sp.Rational(7,10)]])
A = (B+C)/2
SIGNS = list(itertools.product((-1,1),repeat=4))
LOWER = sp.Rational(0)
UPPER = sp.Rational(2,13)
OLD = (-sp.Integer(10413401)+7*sp.sqrt(597925520449617))/1795734335
assert OLD > sp.Rational(1,12)

def exact_certificate(radius, P):
    """Return (bool, diagnostic), with exact rational comparisons only."""
    if not (P == P.T and P[0,0]>0 and sp.factor(P.det())>0):
        return False, "P not positive definite"
    for signs in SIGNS:
        M=A+radius*sp.Matrix(2,2,signs)
        H=P-M.T*P*M
        h11=sp.factor(H[0,0]); h22=sp.factor(H[1,1])
        determinant=sp.factor(H.det())
        if not (h11>0 and h22>0 and determinant>0):
            return False, (signs,h11,h22,determinant)
    return True, "all 16 strict rational vertex inequalities passed"

def sdp_at(radius):
    """Numerical candidate only; solver status is not a proof."""
    matrices=[np.array(A+radius*sp.Matrix(2,2,signs),dtype=float)
              for signs in SIGNS]
    P=cp.Variable((2,2),symmetric=True)
    # trace normalization and modest strict slack improve conditioning.
    constraints=[cp.trace(P)==1, P >> 1e-7*np.eye(2)]
    constraints += [P-M.T@P@M >> 1e-8*np.eye(2) for M in matrices]
    problem=cp.Problem(cp.Minimize(0),constraints)
    try:
        problem.solve(solver="CLARABEL",verbose=False)
    except Exception as exc:
        return None, "solver error: "+str(exc)
    if problem.status not in ("optimal","optimal_inaccurate") or P.value is None:
        return None, problem.status
    return P.value, problem.status

print("Checkpoint 050: numerical SDP exploration; exact certificate required")
lo=float(OLD)
hi=float(UPPER)
best=None
for iteration in range(18):
    mid=(lo+hi)/2
    candidate,status=sdp_at(sp.Rational(str(mid)))
    print("iteration",iteration+1,"radius",format(mid,".12f"),"status",status)
    if candidate is None:
        hi=mid
    else:
        lo=mid
        best=(mid,candidate)
print("NUMERICAL BRACKET ONLY:",lo,hi)

# Rational reconstruction at conservative radii; a failed reconstruction
# is inconclusive, not evidence that no certificate exists.
if best is None:
    print("NO NUMERICAL CANDIDATE; no new mathematical claim")
    sys.exit(0)
for denom in (1000,10000,100000,1000000,10000000):
    # Test at rational radius below numerical candidate, preserving slack.
    target=sp.Rational(str(max(float(OLD)+1e-7,best[0]-0.0005)))
    matrix=best[1]
    p=sp.Rational(str(float(matrix[0,0]))).limit_denominator(denom)
    q=sp.Rational(str(float((matrix[0,1]+matrix[1,0])/2))).limit_denominator(denom)
    r=sp.Rational(1)-p
    Q=sp.Matrix([[p,q],[q,r]])
    passed,detail=exact_certificate(target,Q)
    if passed:
        print("EXACT CERTIFICATE PASS")
        print("radius =",target)
        print("P =",Q)
        print("checked =",detail)
        if target>OLD:
            print("CERTIFIED IMPROVEMENT OVER CHECKPOINT 049")
        else:
            print("No certified improvement over checkpoint 049")
        break
else:
    print("NO EXACT CERTIFICATE RECONSTRUCTED")
    print("Numerical results remain exploratory; do not claim a new bound.")
