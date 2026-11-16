import sympy as sp

print("Checkpoint 036: Active-vertex spectral bound")

e = sp.symbols("e", real=True)

B = sp.Matrix([
    [sp.Rational(4, 5), -2],
    [0, sp.Rational(3, 10)]
])

D = sp.Matrix([
    [1, -1],
    [-1, 1]
])

M = B + e * D

critical = sp.Rational(7, 145)

# Exact characteristic polynomial at lambda = 1.
expression = sp.factor((sp.eye(2) - M).det())
expected = (7 - 145*e)/50

assert sp.simplify(expression - expected) == 0

print("PASS: exact determinant identity")
print("det(I-M) =", expression)

# Exact eigenvalues at the critical perturbation.
Mcrit = M.subs(e, critical)
eigenvalues = Mcrit.eigenvals()

print("critical_epsilon =", critical)
print("critical_decimal =", sp.N(critical, 16))
print("boundary_eigenvalues =", eigenvalues)

assert sp.Integer(1) in eigenvalues
assert max(abs(complex(sp.N(x))) for x in eigenvalues) == 1

print("PASS: unit eigenvalue at boundary")

# Check the spectral transition using exact arithmetic.
for eps in (
    sp.Rational(48, 1000),
    critical,
    sp.Rational(49, 1000)
):
    A = M.subs(e, eps)
    rho = max(abs(complex(sp.N(x))) for x in A.eigenvals())

    print(
        f"epsilon={float(eps):.12f} "
        f"spectral_radius={rho:.12f}"
    )

    if eps < critical:
        assert rho < 1
    elif eps > critical:
        assert rho > 1
    else:
        assert abs(rho - 1) < 1e-12

print("PASS: spectral transition")

# Compare against established quadratic certificate.
certified = sp.Rational(4827, 100000)

assert certified < critical

print("certified_lower =", sp.N(certified, 16))
print("spectral_upper =", sp.N(critical, 16))
print("remaining_gap =", sp.N(critical-certified, 16))

print("PASS: Checkpoint 036 verified")
