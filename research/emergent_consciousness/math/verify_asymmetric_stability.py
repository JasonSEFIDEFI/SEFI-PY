from math import exp

X0 = (1.0, -1.0)
T = 1.0

def rhs(x, a, b):
    d = x[0] - x[1]
    return (-a*d, b*d)

def rk4(x, a, b, dt):
    def add(u, v, c):
        return tuple(i+c*j for i,j in zip(u,v))

    k1 = rhs(x,a,b)
    k2 = rhs(add(x,k1,dt/2),a,b)
    k3 = rhs(add(x,k2,dt/2),a,b)
    k4 = rhs(add(x,k3,dt),a,b)

    return tuple(
        x[i]+dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6
        for i in range(2)
    )

def exact(t,a,b):
    s = a+b
    d0 = X0[0]-X0[1]

    if s == 0:
        factor = t
    else:
        factor = -__import__("math").expm1(-s*t)/s

    return (
        X0[0]-a*d0*factor,
        X0[1]+b*d0*factor
    )

def verify(a,b,n=1000):
    x = X0
    dt = T/n
    conserved = b*X0[0]+a*X0[1]
    max_error = 0.0
    max_conservation_error = 0.0

    for step in range(n+1):
        t = step*dt
        target = exact(t,a,b)

        max_error = max(
            max_error,
            *(abs(x[i]-target[i]) for i in range(2))
        )

        max_conservation_error = max(
            max_conservation_error,
            abs(b*x[0]+a*x[1]-conserved)
        )

        if step < n:
            x = rk4(x,a,b,dt)

    s = a+b
    initial_difference = X0[0]-X0[1]
    final_difference = x[0]-x[1]
    predicted_difference = initial_difference*exp(-s*T)

    assert max_error < 1e-8
    assert max_conservation_error < 1e-8
    assert abs(final_difference-predicted_difference) < 1e-8

    # Check eigenvalue identities via trace and determinant.
    trace = -a-b
    determinant = (-a)*(-b)-a*b
    assert trace == -s
    assert determinant == 0

    if s > 0:
        assert abs(final_difference) < abs(initial_difference)
    elif s < 0:
        assert abs(final_difference) > abs(initial_difference)
    else:
        assert abs(final_difference-initial_difference) < 1e-8

    return x, max_error, max_conservation_error

print("Checkpoint 012: Asymmetric coupling")

cases = [
    (2.0, 1.0),
    (1.0, 3.0),
    (2.0, -1.0),
    (1.0, -1.0),
    (-2.0, 1.0),
]

for a,b in cases:
    x,error,conservation = verify(a,b)

    print(
        f"a={a:+.1f} b={b:+.1f} "
        f"sum={a+b:+.1f} "
        f"x=({x[0]:+.6f},{x[1]:+.6f}) "
        f"error={error:.3e} "
        f"conservation={conservation:.3e}"
    )

print("PASS: asymmetric stability verified")
