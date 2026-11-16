import sympy as sp
from math import exp, sqrt

print("Checkpoint 014: Non-normal transient growth")

A = sp.Matrix([[-1, 10], [0, -2]])
x0 = (1.0, 1.0)

# Exact spectral and non-normality checks
assert A.eigenvals() == {-1: 1, -2: 1}
assert A*A.T != A.T*A

# Symmetric part determines instantaneous
# Euclidean norm growth.
S = (A + A.T)/2
assert S == sp.Matrix([[-1, 5], [5, -2]])
assert S.det() < 0

v = sp.Matrix([1, 1])
initial_growth = 2*(v.T*S*v)[0]
assert initial_growth == 14

print("PASS: exact spectral and growth identities")

def exact(t):
    return (
        11*exp(-t)-10*exp(-2*t),
        exp(-2*t)
    )

def norm(x):
    return sqrt(sum(v*v for v in x))

def rhs(x):
    return (-x[0]+10*x[1], -2*x[1])

def rk4(x, dt):
    def add(a,b,s):
        return tuple(u+s*v for u,v in zip(a,b))

    a = rhs(x)
    b = rhs(add(x,a,dt/2))
    c = rhs(add(x,b,dt/2))
    d = rhs(add(x,c,dt))

    return tuple(
        x[i]+dt*(a[i]+2*b[i]+2*c[i]+d[i])/6
        for i in range(2)
    )

dt = 0.001
steps = 10000
x = x0
max_error = 0.0
max_norm = norm(x)
time_of_max = 0.0

for step in range(steps+1):
    t = step*dt
    target = exact(t)

    max_error = max(
        max_error,
        *(abs(x[i]-target[i]) for i in range(2))
    )

    current_norm = norm(x)
    if current_norm > max_norm:
        max_norm = current_norm
        time_of_max = t

    if step < steps:
        x = rk4(x,dt)

assert max_error < 1e-8
assert max_norm > norm(x0)
assert norm(x) < 0.01

print(f"initial_norm={norm(x0):.6f}")
print(f"maximum_norm={max_norm:.6f}")
print(f"time_of_maximum={time_of_max:.3f}")
print(f"final_norm={norm(x):.6f}")
print(f"max_error={max_error:.3e}")

print("PASS: non-normal transient growth verified")
