import math
import numpy as np

print("Checkpoint 019: Irregular switching")

def cycle_matrix(tau):
    s = 10.0*tau
    return math.exp(-2*tau)*np.array([
        [1.0, s],
        [s, 1.0+s*s]
    ])

def multiplier(tau):
    return math.exp(
        2*(math.asinh(5*tau)-tau)
    )

def log_derivative(tau):
    return 10/math.sqrt(1+25*tau*tau)-2

# Verify exact operator-norm identity.
for tau in (0.1, 1.0, 3.0, 4.0, 5.0, 8.0):
    M = cycle_matrix(tau)
    norm = np.linalg.norm(M, 2)
    expected = multiplier(tau)

    assert abs(norm-expected) < 1e-10

print("PASS: exact operator-norm identity")

# Analytical sufficient condition:
# for tau >= 4,
# log(mu)' <= 2/tau - 2 <= -1.5.
assert log_derivative(4.0) < 0
assert multiplier(4.0) < 1

q = multiplier(4.0)

print(f"uniform_contraction_bound={q:.9f}")
print("PASS: uniform contraction condition")

# Nonperiodic sequence of cycle durations.
schedule = [
    4.0, 5.5, 4.2, 7.0,
    4.8, 6.1, 4.1, 8.0,
    5.0, 4.4, 6.7, 4.3
]

x = np.array([1.0, 1.0])
initial_norm = np.linalg.norm(x)
previous_norm = initial_norm

for n, tau in enumerate(schedule, 1):
    assert tau >= 4.0

    x = cycle_matrix(tau) @ x
    current_norm = np.linalg.norm(x)

    # Bound for this specific cycle.
    assert current_norm <= (
        multiplier(tau)*previous_norm + 1e-12
    )

    # Uniform bound for every allowed schedule.
    assert current_norm <= (
        q**n*initial_norm + 1e-12
    )

    previous_norm = current_norm

    print(
        f"cycle={n:2d} "
        f"tau={tau:.1f} "
        f"norm={current_norm:.9e} "
        f"bound={q**n*initial_norm:.9e}"
    )

print("PASS: irregular switching contraction")

# Contrast: repeated short cycles.
short_M = cycle_matrix(0.1)
eigenvalues, eigenvectors = np.linalg.eigh(short_M)
dominant = eigenvectors[:, -1]

assert eigenvalues[-1] > 1

x = dominant.copy()

for _ in range(8):
    x = short_M @ x

expected_growth = multiplier(0.1)**8

assert np.isclose(
    np.linalg.norm(x),
    expected_growth,
    rtol=1e-11
)

print(f"short_cycle_growth={np.linalg.norm(x):.6f}")
print("PASS: unstable short-cycle contrast")
print("PASS: Checkpoint 019 verified")
