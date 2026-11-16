import math
import sympy as sp
import numpy as np

print("Checkpoint 020: Unequal phase durations")

u, v = sp.symbols("u v", positive=True)

A1 = sp.Matrix([[-1,10],[0,-1]])
A2 = sp.Matrix([[-1,0],[10,-1]])

M = (A2*v).exp() * (A1*u).exp()

expected = sp.exp(-(u+v))*sp.Matrix([
    [1,10*u],
    [10*v,1+100*u*v]
])

assert sp.simplify(M-expected) == sp.zeros(2)
assert sp.simplify(
    M.det()-sp.exp(-2*(u+v))
) == 0

print("PASS: exact unequal-phase operator")

def matrix(a,b):
    factor = math.exp(-(a+b))
    return factor*np.array([
        [1.0,10*a],
        [10*b,1.0+100*a*b]
    ])

def log_multiplier(a,b):
    return -(a+b)+2*math.asinh(
        5*math.sqrt(a*b)
    )

# Independently compare the analytic formula
# against numerical matrix eigenvalues.
durations = [0.1,0.5,1.0,2.0,4.0,6.0]
tested = 0

for a in durations:
    for b in durations:
        eig = np.linalg.eigvals(matrix(a,b))
        radius = max(abs(eig))
        predicted = math.exp(log_multiplier(a,b))

        assert np.isclose(
            radius,predicted,
            rtol=1e-11,atol=1e-12
        )

        # Independent stability condition.
        boundary = (
            5*math.sqrt(a*b)
            -math.sinh((a+b)/2)
        )

        if abs(boundary) > 1e-10:
            assert (radius > 1) == (boundary > 0)

        tested += 1

print("matrix_pairs_tested =",tested)
print("PASS: eigenvalue and boundary agreement")

# Unequal-duration examples.
examples = [
    (0.1,1.0),
    (0.1,8.0),
    (1.0,3.0),
    (3.0,4.0),
    (4.0,5.0)
]

for a,b in examples:
    log_rho = log_multiplier(a,b)
    rho = math.exp(log_rho)
    status = "UNSTABLE" if log_rho > 0 else "STABLE"

    print(
        f"u={a:.2f} v={b:.2f} "
        f"multiplier={rho:.8f} "
        f"{status}"
    )

# Equal-duration consistency with
# Checkpoint 018.
critical = 3.578234808447358

assert abs(
    log_multiplier(critical,critical)
) < 1e-10

print("PASS: equal-phase threshold recovered")
print("PASS: Checkpoint 020 verified")
