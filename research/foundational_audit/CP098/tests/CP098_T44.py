"""CP098-T44: Hamiltonian-flow consistency audit.

Assembles the independently computed T41 bracket matrix
with the exact off-surface Hamiltonian identity.

Checks the T34 normal multiplier independently of
the original preservation-equation assembly.

This is a weak constraint-surface calculation.
Gauge generator analysis remains open.
"""

import sympy as sp

v, k = sp.symbols("v k", positive=True)
A, C = sp.symbols("A C", real=True, nonzero=True)
a, n, r, s, x, y, DA = sp.symbols(
    "a n r s x y DA", real=True
)

LT = sp.symbols("lambda_T", real=True)

d = 2*A/v**3
t = k*v

R = r*r + s*s

# Independently obtained T41 bracket:
E = -(
    4*A*C*a*k
    - 4*C*DA*k*v**2
    + 2*C*n*v**2
    + a*k*R*v**4
    + k*v**5*(r*x+s*y)
)/(2*C*v**3)

# T34 normal preservation forcing:
FN = (
    -n
    - 6*A*a*k/v**2
    - v**3*k*(r*x+s*y)/(2*C)
    - v**2*k*a*R/(2*C)
    + 2*k*DA
)

names = [
    "phi_T",
    "phi_N",
    "chi_T",
    "chi_N",
    "psi_T"
]

M = sp.Matrix([
    [0, 0, 0, 0, 0],
    [0, 0, 0, -d, t*d],
    [0, 0, 0, 2*t*d, -2*t*t*d],
    [0, d, -2*t*d, 0, E],
    [0, -t*d, 2*t*t*d, -E, 0]
])

print("=== CP098-T44 ===", flush=True)

assert M + M.T == sp.zeros(5)
print("PASS: antisymmetric constraint matrix")

# Hc = v*psi_T + v^2*k*chi_N + a*chi_T.
#
# The coefficient derivatives multiply constraints and
# vanish when restricted to the constraint surface.

Hcoeff = sp.Matrix([
    0,
    0,
    a,
    v**2*k,
    v
])

flow = M*Hcoeff

print("Weak canonical Hamiltonian flow:", flush=True)

for name, value in zip(names, flow):
    print(
        "{" + name + ",Hc} =",
        sp.factor(value),
        flush=True
    )

# Compare the chi_N forcing from T34.
forcing = sp.factor(flow[3])

print("Matrix-derived F_N =", forcing, flush=True)

residual = sp.factor(forcing-FN)

print("T34 forcing residual =", residual, flush=True)
assert residual == 0

print("PASS: independent F_N agreement")

# Primary multipliers in the original basis:
# Htotal = Hc + LT*phi_T + LN*phi_N

LN = -FN/d

multipliers = sp.Matrix([
    LT,
    LN,
    0,
    0,
    0
])

total_flow = sp.simplify(flow+M*multipliers)

print("Normal multiplier =", sp.factor(LN),
      flush=True)

for name, value in zip(names, total_flow):
    residual = sp.factor(value)
    print(
        "dot(" + name + ") =",
        residual,
        flush=True
    )
    assert residual == 0, name

print("PASS: all five total preservation equations")

# LT should disappear identically.
assert not any(
    LT in expression.free_symbols
    for expression in total_flow
)

print("PASS: lambda_T remains unrestricted")
print("PASS: T34 and T41 agree on lambda_N")
print("REGULAR HAMILTONIAN FLOW: CONSISTENT")
print("GAUGE GENERATOR: NOT YET DERIVED")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
