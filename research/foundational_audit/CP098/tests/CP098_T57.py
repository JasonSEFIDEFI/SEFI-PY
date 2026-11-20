"""CP098-T57: Derive first-class normal corrections.

Solve the weak first-class conditions for a general
five-constraint generator.

The coefficient functions may depend on phase space.
Their Poisson derivatives multiply constraints and
therefore vanish on the regular constraint surface.

No action invariance or gauge generator is assumed.
"""

import sympy as sp

v, k = sp.symbols("v k", positive=True)
A, C = sp.symbols("A C", real=True, nonzero=True)

a, n, r, s, x, y, DA = sp.symbols(
    "a n r s x y DA", real=True
)

f, g, h, j, ell = sp.symbols(
    "f g h j ell", real=True
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

# Constraint order:
# phi_T, phi_N, chi_T, chi_N, psi_T

M = sp.Matrix([
    [0, 0, 0, 0, 0],
    [0, 0, 0, -d, t*d],
    [0, 0, 0, 2*t*d, -2*t*t*d],
    [0, d, -2*t*d, 0, E],
    [0, -t*d, 2*t*t*d, -E, 0]
])

print("=== CP098-T57 ===", flush=True)

assert M + M.T == sp.zeros(5)
assert M.rank() == 2

print("PASS: regular rank-two matrix", flush=True)

# General generator coefficient vector.
coeff = sp.Matrix([f, j, g, ell, h])

# {G,constraint_i} weakly equals
# sum_j coeff_j * {constraint_j,constraint_i}.
#
# This is -M*coeff because M is antisymmetric.

residual = (-M*coeff).applyfunc(sp.factor)

names = [
    "phi_T",
    "phi_N",
    "chi_T",
    "chi_N",
    "psi_T"
]

print("Weak first-class residuals:")
for name, value in zip(names, residual):
    print(name, ":", value, flush=True)

# Solve the two independent conditions.
#
# Against phi_N:
# d*(ell-t*h) = 0.
#
# Against chi_N:
# -d*j + 2*t*d*g - E*h = 0.

normal_primary_residual = sp.factor(residual[1])
normal_secondary_residual = sp.factor(residual[3])

assert sp.factor(
    normal_primary_residual-d*(ell-t*h)
) == 0

assert sp.factor(
    normal_secondary_residual
    -(-d*j+2*t*d*g-E*h)
) == 0

print("PASS: independent normal conditions",
      flush=True)

solution = sp.solve(
    [
        normal_primary_residual,
        normal_secondary_residual
    ],
    [ell, j],
    dict=True
)

assert len(solution) == 1

derived = solution[0]

print("Derived ell =", sp.factor(derived[ell]),
      flush=True)
print("Derived j =", sp.factor(derived[j]),
      flush=True)

assert sp.factor(derived[ell]-t*h) == 0
assert sp.factor(
    derived[j]-(2*t*g-E*h/d)
) == 0

print("PASS: normal corrections uniquely fixed",
      flush=True)

# Substitute into every first-class condition.
corrected = residual.subs(derived)
corrected = corrected.applyfunc(sp.factor)

assert corrected == sp.zeros(5, 1)

print("PASS: all five weak brackets vanish",
      flush=True)

# Reconstruct the previously identified basis:
#
# G1 = phi_T
# G2 = chi_T + 2*t*phi_N
# G3 = psi_T + t*chi_N - (E/d)*phi_N

G1 = sp.Matrix([1, 0, 0, 0, 0])
G2 = sp.Matrix([0, 2*t, 1, 0, 0])
G3 = sp.Matrix([0, -E/d, 0, t, 1])

reconstruction = f*G1+g*G2+h*G3

corrected_coeff = coeff.subs(derived)

assert all(
    sp.factor(corrected_coeff[i]-reconstruction[i]) == 0
    for i in range(5)
)

print("PASS: T42 basis independently recovered",
      flush=True)

# Verify that psi_T by itself is not first-class.
bare_tertiary = sp.Matrix([0, 0, 0, 0, 1])
bare_residual = (-M*bare_tertiary).applyfunc(
    sp.factor
)

assert bare_residual != sp.zeros(5, 1)
assert sp.factor(bare_residual[1]+t*d) == 0

print("PASS: bare tertiary is not first-class",
      flush=True)

# Verify that the corrected tertiary is first-class.
assert M*G3 == sp.zeros(5, 1)

print("PASS: corrected tertiary first-class weakly",
      flush=True)

# Rank of the reconstructed basis.
basis = sp.Matrix.hstack(G1, G2, G3)

assert basis.rank() == 3

print("PASS: three independent first-class directions",
      flush=True)

# The original primary first-class intersection
# remains one-dimensional.
assert M*G1 == sp.zeros(5, 1)

normal_primary = sp.Matrix([0, 1, 0, 0, 0])
assert M*normal_primary != sp.zeros(5, 1)

print("PASS: primary first-class dimension remains 1",
      flush=True)

print("GENERAL WEAK FIRST-CLASS GENERATOR:")
print("G = f*G1 + g*G2 + h*G3")

print("NORMAL CORRECTIONS:")
print("ell = t*h")
print("j = 2*t*g - (E/d)*h")

print("WEAK FIRST-CLASS CORRECTIONS: VERIFIED")
print("TIME-DEPENDENT COEFFICIENT RELATIONS: OPEN")
print("OFF-SHELL ACTION INVARIANCE: OPEN")
print("MULTIPLIER TRANSFORMATIONS: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
