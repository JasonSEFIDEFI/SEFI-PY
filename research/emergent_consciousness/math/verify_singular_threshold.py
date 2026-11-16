import math
import numpy as np
import sympy as sp

print("Checkpoint 023: Singular-value threshold")

s = 0.01

# Exact symbolic identities
a, L = sp.symbols("a L", positive=True)
Q = sp.Matrix([
    [1, 10*a],
    [10*L, 1+100*a*L]
])

F = sp.expand(sp.trace(Q.T*Q))

assert sp.expand(Q.det()-1) == 0
assert sp.expand(
    F-(2+100*a*a+100*L*L+200*a*L
       +10000*a*a*L*L)
) == 0

print("PASS: exact singular-value identities")

def matrix(length):
    return math.exp(-(s+length))*np.array([
        [1.0,10*s],
        [10*length,1+100*s*length]
    ])

def boundary(length):
    f = (2+100*s*s+100*length*length
         +200*s*length
         +10000*s*s*length*length)
    return f-2*math.cosh(2*(s+length))

# Locate the root in the established bracket.
left, right = 3.0, 4.0
assert boundary(left) > 0
assert boundary(right) < 0

for _ in range(80):
    middle = (left+right)/2
    if boundary(middle) > 0:
        left = middle
    else:
        right = middle

critical = (left+right)/2

assert abs(boundary(critical)) < 1e-9
assert abs(critical-3.575084101323) < 1e-10

print(f"critical_L={critical:.12f}")
print("PASS: threshold root located")

# Compare independent SVD against the
# analytical squared-singular-value formula.
for length in (
    3.0, 3.5, critical-0.01,
    critical, critical+0.01, 4.0
):
    B = matrix(length)
    singular = np.linalg.svd(
        B, compute_uv=False
    )[0]

    f = (2+100*s*s+100*length*length
         +200*s*length
         +10000*s*s*length*length)

    predicted_squared = math.exp(
        -2*(s+length)
    )*(f+math.sqrt(f*f-4))/2

    assert np.isclose(
        singular**2,
        predicted_squared,
        rtol=1e-12, atol=1e-12
    )

    print(
        f"L={length:.9f} "
        f"sigma_max={singular:.12f}"
    )

print("PASS: independent SVD agreement")

# Verify classification on both sides.
for length in (critical-0.01, critical+0.01):
    B = matrix(length)
    C = B.T
    sigma = np.linalg.norm(B, 2)
    product_radius = max(
        abs(np.linalg.eigvals(C@B))
    )

    assert np.isclose(
        product_radius, sigma**2,
        rtol=1e-12
    )

    if length < critical:
        assert sigma > 1
        assert product_radius > 1
    else:
        assert sigma < 1
        assert product_radius < 1

print("PASS: boundary classification")
print("PASS: Checkpoint 023 verified")
