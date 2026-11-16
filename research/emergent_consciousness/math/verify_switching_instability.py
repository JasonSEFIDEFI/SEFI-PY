import sympy as sp
from math import exp, sqrt

print("Checkpoint 017: Switching instability")

A1 = sp.Matrix([[-1,10],[0,-1]])
A2 = sp.Matrix([[-1,0],[10,-1]])

tau = sp.Rational(1,10)

assert A1.eigenvals() == {-1:2}
assert A2.eigenvals() == {-1:2}

E1 = (A1*tau).exp()
E2 = (A2*tau).exp()
M = sp.simplify(E2*E1)

expected = sp.exp(-sp.Rational(1,5))*sp.Matrix([
    [1,1],
    [1,2]
])

assert sp.simplify(M-expected) == sp.zeros(2)

mu_plus = sp.exp(-sp.Rational(1,5))*(
    3+sp.sqrt(5)
)/2

mu_minus = sp.exp(-sp.Rational(1,5))*(
    3-sp.sqrt(5)
)/2

assert abs(float(mu_plus)-2.143) < 0.002
assert float(mu_plus) > 1
assert float(mu_minus) < 1

print("PASS: exact switching operator")
print(f"cycle_eigenvalues=({float(mu_plus):.6f},"
      f" {float(mu_minus):.6f})")

def rhs(x, matrix):
    return (
        matrix[0][0]*x[0]+matrix[0][1]*x[1],
        matrix[1][0]*x[0]+matrix[1][1]*x[1]
    )

def rk4(x, matrix, dt):
    def add(u,v,s):
        return tuple(a+s*b for a,b in zip(u,v))

    a = rhs(x,matrix)
    b = rhs(add(x,a,dt/2),matrix)
    c = rhs(add(x,b,dt/2),matrix)
    d = rhs(add(x,c,dt),matrix)

    return tuple(
        x[i]+dt*(a[i]+2*b[i]+2*c[i]+d[i])/6
        for i in range(2)
    )

def matvec(matrix,x):
    return (
        matrix[0][0]*x[0]+matrix[0][1]*x[1],
        matrix[1][0]*x[0]+matrix[1][1]*x[1]
    )

def norm(x):
    return sqrt(x[0]**2+x[1]**2)

B1 = ((-1.0,10.0),(0.0,-1.0))
B2 = ((-1.0,0.0),(10.0,-1.0))

factor = exp(-0.2)
cycle = (
    (factor,factor),
    (factor,2*factor)
)

x0 = (1.0,1.0)
numeric = x0
reference = x0

dt = 0.001
steps_per_phase = 100
cycles = 12
max_relative_error = 0.0

for n in range(cycles):
    for matrix in (B1,B2):
        for _ in range(steps_per_phase):
            numeric = rk4(numeric,matrix,dt)

    reference = matvec(cycle,reference)

    relative_error = norm((
        numeric[0]-reference[0],
        numeric[1]-reference[1]
    ))/max(1.0,norm(reference))

    max_relative_error = max(
        max_relative_error,relative_error
    )

    assert relative_error < 1e-8

    if n in (0,1,3,7,11):
        print(
            f"cycle={n+1:2d} "
            f"norm={norm(numeric):.6f} "
            f"relative_error={relative_error:.3e}"
        )

assert norm(numeric) > 1000*norm(x0)
assert max_relative_error < 1e-8

print(f"initial_norm={norm(x0):.6f}")
print(f"final_norm={norm(numeric):.6f}")
print(f"max_relative_error={max_relative_error:.3e}")
print("PASS: switching instability verified")
