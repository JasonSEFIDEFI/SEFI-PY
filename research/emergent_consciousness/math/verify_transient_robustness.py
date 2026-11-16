import sympy as sp
from math import sqrt

print("Checkpoint 015: Structural robustness")

c = sp.symbols("c", real=True)

A = sp.Matrix([[-1, c], [0, -2]])
S = (A + A.T)/2

# Characteristic polynomial of the symmetric part.
z = sp.symbols("z")
expected = z*z + 3*z + 2 - c*c/4

assert sp.expand(
    S.charpoly(z).as_expr().subs(
        S.charpoly(z).gen, z
    ) - expected
) == 0

print("PASS: exact symmetric-part polynomial")

def growth_rate(coupling, x):
    a,b = x
    return 2*(-a*a + coupling*a*b - 2*b*b)

def largest_symmetric_eigenvalue(coupling):
    return (-3 + sqrt(1+coupling*coupling))/2

threshold = 2*sqrt(2)

for coupling in (0.0, 2.0, threshold, 3.0, 5.0, 10.0):
    eigenvalue = largest_symmetric_eigenvalue(coupling)

    if coupling < threshold - 1e-10:
        assert eigenvalue < 0
    elif coupling > threshold + 1e-10:
        assert eigenvalue > 0
    else:
        assert abs(eigenvalue) < 1e-10

    # Search a dense set of unit-vector directions.
    from math import cos, sin, pi
    maximum_growth = max(
        growth_rate(
            coupling,
            (cos(2*pi*j/10000), sin(2*pi*j/10000))
        )
        for j in range(10000)
    )

    # Maximum rate for a unit vector is 2*lambda_max(S).
    assert abs(maximum_growth-2*eigenvalue) < 1e-6

    print(
        f"c={coupling:.6f}, "
        f"lambda_max={eigenvalue:+.6f}, "
        f"max_growth={maximum_growth:+.6f}"
    )

print("PASS: transient-growth threshold verified")
