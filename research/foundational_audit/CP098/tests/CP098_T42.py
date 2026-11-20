"""CP098-T42: Algebraic first-class basis audit.

Uses the T41 bracket-matrix structure.
Checks three independent null combinations and
a nondegenerate second-class pair.

This is a matrix-level verification, not an
independent recomputation of T41's brackets.

Gauge generator and singular sectors remain open.
"""

import sympy as sp

d, t, E = sp.symbols(
    "d t E", real=True, nonzero=True
)

# Constraint ordering:
# phi_T, phi_N, chi_T, chi_N, psi_T

M = sp.Matrix([
    [0, 0, 0, 0, 0],
    [0, 0, 0, -d, t*d],
    [0, 0, 0, 2*t*d, -2*t**2*d],
    [0, d, -2*t*d, 0, E],
    [0, -t*d, 2*t**2*d, -E, 0],
])

G1 = sp.Matrix([1, 0, 0, 0, 0])

# chi_T + 2*t*phi_N
G2 = sp.Matrix([0, 2*t, 1, 0, 0])

# psi_T + t*chi_N - (E/d)*phi_N
G3 = sp.Matrix([0, -E/d, 0, t, 1])

G = sp.Matrix.hstack(G1, G2, G3)

print("=== CP098-T42 ===", flush=True)

assert M + M.T == sp.zeros(5)
print("PASS: antisymmetry", flush=True)

for label, vec in [
    ("G1", G1),
    ("G2", G2),
    ("G3", G3),
]:
    residual = sp.simplify(M*vec)
    print(label, "null residual =", list(residual),
          flush=True)
    assert residual == sp.zeros(5, 1)

print("PASS: three null combinations", flush=True)

# Independence follows from the minor formed
# by rows 0, 2, and 4.
independence_minor = sp.factor(
    G.extract([0, 2, 4], [0, 1, 2]).det()
)

print("Independence minor =", independence_minor,
      flush=True)
assert independence_minor != 0

print("PASS: independent null directions",
      flush=True)

# phi_N and chi_N form a second-class pair.
second_class = M.extract([1, 3], [1, 3])
det_sc = sp.factor(second_class.det())

print("Second-class submatrix =", second_class,
      flush=True)
print("Second-class determinant =", det_sc,
      flush=True)

assert sp.simplify(det_sc-d**2) == 0

# Rank <= 2: all columns are spanned by
# columns 1 and 3.
basis = M[:, [1, 3]]
for j in range(5):
    # Direct rank check on compact matrix.
    assert sp.Matrix.hstack(basis, M[:, j]).rank() <= 2

assert M.rank() == 2

print("PASS: exact matrix rank = 2", flush=True)
print("FIRST-CLASS COMBINATIONS: 3")
print("SECOND-CLASS CONSTRAINTS: 2")
print("GAUGE GENERATOR: NOT YET DERIVED")
print("SINGULAR-SECTOR CLASSIFICATION: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
