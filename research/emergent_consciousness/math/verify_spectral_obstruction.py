import sympy as sp
import numpy as np

print("Checkpoint 031: Spectral obstruction")

e = sp.symbols("e", real=True)

B = sp.Matrix([
    [sp.Rational(4,5), -2],
    [0, sp.Rational(3,10)]
])

E = e*sp.Matrix([
    [1,1],
    [-1,-1]
])

M = B+E

assert sp.simplify(
    sp.trace(M)-sp.Rational(11,10)
) == 0

assert sp.simplify(
    M.det()-(sp.Rational(6,25)
             -sp.Rational(5,2)*e)
) == 0

characteristic_at_one = sp.factor(
    (sp.eye(2)-M).det()
)

expected = (
    sp.Rational(7,50)
    -sp.Rational(5,2)*e
)

assert sp.simplify(
    characteristic_at_one-expected
) == 0

critical = sp.Rational(7,125)

assert sp.simplify(
    characteristic_at_one.subs(e,critical)
) == 0

print("PASS: exact spectral boundary")
print("critical_epsilon =", critical)
print("decimal =", float(critical))

# At the boundary, verify eigenvalues exactly.
Mcrit = M.subs(e,critical)
eigenvalues = Mcrit.eigenvals()

assert sp.Integer(1) in eigenvalues
assert sp.Rational(1,10) in eigenvalues

print("PASS: exact unit eigenvalue")
print("boundary_eigenvalues =", list(eigenvalues))

# Numerical checks on either side.
for value in (0.05,0.056,0.06):
    A = np.array(M.subs(e,value),dtype=float)
    rho = max(abs(np.linalg.eigvals(A)))

    print(
        f"epsilon={value:.6f} "
        f"spectral_radius={rho:.12f}"
    )

    if value < float(critical):
        assert rho < 1
    elif value > float(critical):
        assert rho > 1
    else:
        assert np.isclose(rho,1,atol=1e-12)

print("PASS: spectral transition verified")

previous = sp.Rational(13,375)

assert previous < critical

print("previous_metric_boundary =", float(previous))
print("spectral_upper_bound =", float(critical))
print("unresolved_gap =", float(critical-previous))

print("PASS: Checkpoint 031 verified")
