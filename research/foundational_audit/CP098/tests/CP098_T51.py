"""CP098-T51: First-class Hamiltonian-flow consistency.

Checks preservation equations in the T42 first-class /
second-class basis, keeping primary multipliers free.

Weak constraint-surface calculation.
Not an independent canonical differentiation.
Not a gauge-generator derivation.
"""

import sympy as sp

v, k = sp.symbols("v k", positive=True)
A, C = sp.symbols("A C", nonzero=True, real=True)

a, n, r, s, x, y, DA = sp.symbols(
    "a n r s x y DA", real=True
)

LT, LN = sp.symbols(
    "lambda_T lambda_N", real=True
)

d = 2*A/v**3
t = v*k
R = r*r+s*s

E = -(
    4*A*C*a*k
    - 4*C*DA*k*v**2
    + 2*C*n*v**2
    + a*k*R*v**4
    + k*v**5*(r*x+s*y)
)/(2*C*v**3)

FN = (
    -n
    - 6*A*a*k/v**2
    - v**3*k*(r*x+s*y)/(2*C)
    - v**2*k*a*R/(2*C)
    + 2*k*DA
)

# T41 matrix, original constraint order:
# phi_T, phi_N, chi_T, chi_N, psi_T

M = sp.Matrix([
    [0, 0, 0, 0, 0],
    [0, 0, 0, -d, t*d],
    [0, 0, 0, 2*t*d, -2*t*t*d],
    [0, d, -2*t*d, 0, E],
    [0, -t*d, 2*t*t*d, -E, 0]
])

assert M + M.T == sp.zeros(5)

# Exact Hamiltonian decomposition, evaluated weakly.
h = sp.Matrix([0, 0, a, v**2*k, v])

# Total Hamiltonian adds only primary constraints.
multipliers = sp.Matrix([LT, LN, 0, 0, 0])

flow = (M*(h+multipliers)).applyfunc(sp.factor)

print("=== CP098-T51 ===", flush=True)

print("Original-basis total flow:")
for name, value in zip(
    ["phi_T", "phi_N", "chi_T", "chi_N", "psi_T"],
    flow
):
    print(name, ":", value, flush=True)

assert sp.factor(flow[3]-(FN+d*LN)) == 0
assert sp.factor(flow[4]+t*(FN+d*LN)) == 0

print("PASS: normal forcing reconstructed", flush=True)

# T42 constraint basis:
# G1 = phi_T
# G2 = chi_T + 2*t*phi_N
# G3 = psi_T + t*chi_N - (E/d)*phi_N
# S1 = phi_N
# S2 = chi_N

S = sp.Matrix([
    [1, 0,    0, 0, 0],
    [0, 2*t,  1, 0, 0],
    [0, -E/d, 0, t, 1],
    [0, 1,    0, 0, 0],
    [0, 0,    0, 1, 0]
])

assert sp.factor(S.det()) != 0

new_flow = (S*flow).applyfunc(sp.factor)

labels = ["G1", "G2", "G3", "S1", "S2"]

print("New-basis total flow:")
for name, value in zip(labels, new_flow):
    print(name, ":", value, flush=True)

# The first-class combinations must have vanishing
# weak evolution before fixing LN.
for i in range(3):
    assert sp.factor(new_flow[i]) == 0

print(
    "PASS: all three first-class combinations "
    "preserved for arbitrary primary multipliers",
    flush=True
)

assert sp.factor(new_flow[3]) == 0
assert sp.factor(new_flow[4]-(FN+d*LN)) == 0

print("PASS: only S2 fixes the normal multiplier",
      flush=True)

LN_fixed = -FN/d

fixed_flow = new_flow.subs(LN, LN_fixed)
fixed_flow = fixed_flow.applyfunc(sp.factor)

assert fixed_flow == sp.zeros(5, 1)

print("PASS: full preservation after fixing LN",
      flush=True)

assert not any(
    LT in value.free_symbols
    for value in new_flow
)

print("PASS: tangential multiplier unrestricted",
      flush=True)

print("FIRST-CLASS BASIS FLOW: VERIFIED WEAKLY")
print("NORMAL SECOND-CLASS MULTIPLIER: FIXED")
print("GAUGE GENERATOR: NOT YET DERIVED")
print("OFF-SURFACE FIRST-CLASS ALGEBRA: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
