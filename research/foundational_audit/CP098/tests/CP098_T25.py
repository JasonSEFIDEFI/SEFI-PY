"""CP098-T25: General internal canonical contribution.

Independent symbolic derivation of the internal part
of {chi_N,H_c}, with coordinate-dependent A and C.

This does not test the external canonical sectors
or establish general Dirac-Bergmann closure.
"""

import sympy as sp

I1, I2 = sp.symbols("I1 I2", real=True)
p1, p2 = sp.symbols("p1 p2", real=True)

h1, h2 = sp.symbols("h1 h2", positive=True)
v, k = sp.symbols("v k", positive=True)
P2, m = sp.symbols("P2 m", real=True)

A = sp.Function("A")(I1, I2)
C = sp.Function("C")(I1, I2)

B = (
    m
    + p1**2/(2*h1)
    + p2**2/(2*h2)
)

# Only internal dependence is relevant here.
chi = 2*A*k/v

H = (
    v*(B + A*k**2)
    - v**5*k**2*P2/(4*C)
)

coordinates = (I1, I2)
momenta = (p1, p2)

direct = sum(
    sp.diff(chi, I)*sp.diff(H, p)
    - sp.diff(chi, p)*sp.diff(H, I)
    for I, p in zip(coordinates, momenta)
)

direct = sp.simplify(direct)

D_A = sum(
    sp.diff(A, I)*p/h
    for I, p, h in zip(
        coordinates, momenta, (h1, h2)
    )
)

expected = 2*k*D_A
difference = sp.simplify(direct - expected)

print("=== CP098-T25 ===")
print("Direct internal bracket =", direct)
print("Derived D_A =", D_A)
print("Expected =", expected)
print("Difference =", difference)

assert difference == 0

print("PASS: general two-coordinate internal identity")
print("External canonical sectors: NOT TESTED")
print("GENERAL DIRAC CLOSURE: NOT ESTABLISHED")
