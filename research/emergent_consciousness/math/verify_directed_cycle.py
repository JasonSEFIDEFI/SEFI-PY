import sympy as sp
from math import exp, sqrt

print("Checkpoint 013: Directed-cycle geometry")

k, lam = sp.symbols("k lam", real=True)

P = sp.Matrix([
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 0]
])

A = k*(P-sp.eye(3))

expected = lam*(lam**2 + 3*k*lam + 3*k**2)
characteristic = A.charpoly()
z = characteristic.gen

actual = characteristic.as_expr()
expected = z*(z**2 + 3*k*z + 3*k**2)

assert sp.expand(actual-expected) == 0
assert P**3 == sp.eye(3)
assert A*sp.ones(3,1) == sp.zeros(3,1)
assert sp.ones(1,3)*A == sp.zeros(1,3)

print("PASS: exact spectral identities")

def rhs(x, rate):
    return (
        rate*(x[1]-x[0]),
        rate*(x[2]-x[1]),
        rate*(x[0]-x[2])
    )

def rk4(x, rate, dt):
    def add(a,b,s):
        return tuple(u+s*v for u,v in zip(a,b))

    a = rhs(x,rate)
    b = rhs(add(x,a,dt/2),rate)
    c = rhs(add(x,b,dt/2),rate)
    d = rhs(add(x,c,dt),rate)

    return tuple(
        x[i]+dt*(a[i]+2*b[i]+2*c[i]+d[i])/6
        for i in range(3)
    )

def exact(x0, rate, t):
    z = rate*t

    # Exponential series grouped by n modulo 3.
    coefficients = [0.0,0.0,0.0]
    term = 1.0

    for n in range(120):
        if n:
            term *= z/n
        coefficients[n%3] += term

    p1 = (x0[1],x0[2],x0[0])
    p2 = (x0[2],x0[0],x0[1])

    return tuple(
        exp(-z)*(
            coefficients[0]*x0[i]
            + coefficients[1]*p1[i]
            + coefficients[2]*p2[i]
        )
        for i in range(3)
    )

def disagreement(x):
    mean = sum(x)/3
    return sqrt(sum((v-mean)**2 for v in x))

x0 = (2.0,-1.0,4.0)
rate = 1.0
dt = 0.001
steps = 1000

x = x0
initial_mean = sum(x0)/3
max_error = 0.0
max_mean_error = 0.0
previous_disagreement = disagreement(x)

for step in range(steps+1):
    t = step*dt
    reference = exact(x0,rate,t)

    max_error = max(
        max_error,
        *(abs(x[i]-reference[i]) for i in range(3))
    )

    max_mean_error = max(
        max_mean_error,
        abs(sum(x)/3-initial_mean)
    )

    current_disagreement = disagreement(x)

    assert current_disagreement <= (
        previous_disagreement + 1e-10
    )

    previous_disagreement = current_disagreement

    if step < steps:
        x = rk4(x,rate,dt)

assert max_error < 1e-8
assert max_mean_error < 1e-10
assert disagreement(x) < disagreement(x0)

print(f"final_state={tuple(round(v,6) for v in x)}")
print(f"max_error={max_error:.3e}")
print(f"mean_error={max_mean_error:.3e}")
print(f"initial_disagreement={disagreement(x0):.6f}")
print(f"final_disagreement={disagreement(x):.6f}")
print("PASS: directed-cycle geometry verified")
