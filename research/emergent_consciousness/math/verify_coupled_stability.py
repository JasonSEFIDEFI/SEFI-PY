from math import exp, isclose

X0 = (1.0, -1.0)
DT = 0.001
STEPS = 1000

def rhs(x, k):
    a, b = x
    return (k * (b-a), k * (a-b))

def rk4(x, k, dt):
    def add(y, z, scale):
        return tuple(a + scale*b for a, b in zip(y, z))

    k1 = rhs(x, k)
    k2 = rhs(add(x, k1, dt/2), k)
    k3 = rhs(add(x, k2, dt/2), k)
    k4 = rhs(add(x, k3, dt), k)

    return tuple(
        a + dt*(b + 2*c + 2*d + e)/6
        for a, b, c, d, e in zip(x, k1, k2, k3, k4)
    )

def analytical(t, k):
    mean = sum(X0)/2
    difference = (X0[0]-X0[1])*exp(-2*k*t)
    return (
        mean + difference/2,
        mean - difference/2
    )

def experiment(k):
    x = X0
    max_error = 0.0
    max_mean_error = 0.0

    for step in range(STEPS + 1):
        t = step * DT
        exact = analytical(t, k)

        max_error = max(
            max_error,
            *(abs(a-b) for a, b in zip(x, exact))
        )

        max_mean_error = max(
            max_mean_error,
            abs(sum(x)/2 - sum(X0)/2)
        )

        if step < STEPS:
            x = rk4(x, k, DT)

    initial_difference = abs(X0[0]-X0[1])
    final_difference = abs(x[0]-x[1])

    assert max_error < 1e-5
    assert max_mean_error < 1e-8

    if k > 0:
        assert final_difference < initial_difference
    elif k == 0:
        assert isclose(final_difference, initial_difference)
    else:
        assert final_difference > initial_difference

    return max_error, max_mean_error, final_difference

print("Checkpoint 011: Coupled stability")

for k in (1.0, 0.0, -1.0):
    result = experiment(k)
    repeated = experiment(k)

    assert result == repeated

    error, mean_error, difference = result

    print(
        f"k={k:+.1f} "
        f"final_difference={difference:.6f} "
        f"max_error={error:.3e} "
        f"mean_error={mean_error:.3e}"
    )

print("PASS: stability criteria verified")
