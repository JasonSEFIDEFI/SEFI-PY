"""CP098-T33: Internal tertiary-bracket audit.

Independent canonical differentiation for two internal
pairs, arbitrary differentiable A(I), C(I), and
quadratic internal kinetic energy.

Does not establish full Dirac-Bergmann closure.
"""

import sympy as sp

I1, I2 = sp.symbols("I1 I2", real=True)
p1, p2 = sp.symbols("p1 p2", real=True)

h1, h2 = sp.symbols("h1 h2", positive=True)
k1, k2, m = sp.symbols("k1 k2 m", real=True)

v, kap = sp.symbols("v kap", positive=True)
P2 = sp.symbols("P2", real=True)

A = sp.Function("A")(I1, I2)
C = sp.Function("C")(I1, I2)

B = (
    m
    + p1**2/(2*h1)
    + p2**2/(2*h2)
    + k1*I1**2/2
    + k2*I2**2/2
)

# Terms independent of internal canonical variables
# contribute zero to the internal Poisson bracket.

psi = (
    B - A*kap**2
    - v**4*kap**2*P2/(4*C)
)

chi = 2*A*kap/v

pairs = [(I1, p1), (I2, p2)]

print("=== CP098-T33 ===", flush=True)

direct = sum(
    sp.diff(psi, I)*sp.diff(chi, p)
    - sp.diff(psi, p)*sp.diff(chi, I)
    for I, p in pairs
)

direct = sp.simplify(direct)

DA = (
    sp.diff(A, I1)*p1/h1
    + sp.diff(A, I2)*p2/h2
)

FN_internal = 2*kap*DA
predicted = -FN_internal/v

difference = sp.simplify(direct-predicted)

print("Direct internal bracket =", direct, flush=True)
print("Predicted =", predicted, flush=True)
print("Difference =", difference, flush=True)

assert difference == 0

# Independently check that chi has no internal
# momentum dependence. This is why derivatives
# of C(I) in psi do not enter this bracket.

for p in (p1, p2):
    assert sp.diff(chi, p) == 0

print("PASS: arbitrary differentiable A and C")
print("PASS: internal tertiary identity")
print("FULL REGULAR-SURFACE IDENTITY: pending combined audit")
print("GENERAL DIRAC CLOSURE: NOT ESTABLISHED")
print("DYNAMICAL CONSISTENCY: OPEN")
