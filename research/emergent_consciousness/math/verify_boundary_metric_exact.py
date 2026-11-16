import itertools
import sympy as sp

print("Checkpoint 039: Exact boundary metric")

e = sp.Rational(7, 145)
S = sp.Matrix([[1, 4], [0, 1]])

B = S * sp.diag(sp.Rational(4, 5),
                sp.Rational(3, 10)) * S.inv()
C = S * sp.diag(sp.Rational(2, 5),
                sp.Rational(7, 10)) * S.inv()

D = sp.Matrix([[1, -1], [-1, 1]])
M1 = B + e * D

nullspace = (M1 - sp.eye(2)).nullspace()
assert len(nullspace) == 1
v = nullspace[0]

print("unit_eigenvector =", v.T)

p, q, r = sp.symbols("p q r")
P = sp.Matrix([[p, q], [q, r]])

H1 = P - M1.T * P * M1
equations = list(H1 * v)

# Normalize trace(P)=1.
solutions = sp.linsolve(
    equations + [p + r - 1],
    (p, q, r)
)

solution = next(iter(solutions))
print("boundary_metric_family =", solution)

free = sorted(
    set().union(*(x.free_symbols for x in solution)),
    key=str
)

assert len(free) == 1
parameter = free[0]

vertices = []
for name, base in (("B", B), ("C", C)):
    for signs in itertools.product((-1, 1), repeat=4):
        M = base + e * sp.Matrix(2, 2, signs)
        vertices.append((name, signs, M))

def psd2(A):
    return (
        A[0, 0] >= 0
        and A[1, 1] >= 0
        and A.det() >= 0
    )

def pd2(A):
    return A[0, 0] > 0 and A.det() > 0

found = None

# Rational grid, exact verification at every point.
for k in range(-2000, 2001):
    value = sp.Rational(k, 10000)
    vals = {
        parameter: value
    }
    entries = [sp.factor(x.subs(vals)) for x in solution]
    Q = sp.Matrix([
        [entries[0], entries[1]],
        [entries[1], entries[2]]
    ])

    if not pd2(Q):
        continue

    if all(
        psd2(Q - M.T * Q * M)
        for _, _, M in vertices
    ):
        found = (value, Q)
        break

if found is None:
    print("No rational grid candidate found")
    print("Checkpoint 039 remains unresolved")
else:
    value, Q = found
    print("parameter =", value)
    print("exact_P =")
    sp.pprint(Q)

    assert pd2(Q)

    for _, _, M in vertices:
        assert psd2(Q - M.T * Q * M)

    print("PASS: P positive definite")
    print("PASS: 32 exact boundary inequalities")
    print("PASS: Checkpoint 039 verified")
