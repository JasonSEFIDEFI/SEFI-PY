import math
import numpy as np

print("Checkpoint 018: Switching threshold")

def cycle_matrix(tau):
    s = 10.0*tau
    factor = math.exp(-2.0*tau)
    return factor*np.array([
        [1.0, s],
        [s, 1.0+s*s]
    ])

def log_multiplier(tau):
    return 2.0*(math.asinh(5.0*tau)-tau)

def exact_multiplier(tau):
    return math.exp(log_multiplier(tau))

# Compare two independent formulations.
for tau in (0.01, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0):
    eigenvalues = np.linalg.eigvalsh(cycle_matrix(tau))
    spectral_radius = max(abs(eigenvalues))
    predicted = exact_multiplier(tau)

    assert abs(spectral_radius-predicted) < 1e-10

print("PASS: exact multiplier identity")

# The nonzero root is bracketed away from tau=0.
left = 1.0
right = 5.0

assert log_multiplier(left) > 0
assert log_multiplier(right) < 0

for _ in range(100):
    middle = (left+right)/2.0
    if log_multiplier(middle) > 0:
        left = middle
    else:
        right = middle

critical = (left+right)/2.0

assert abs(log_multiplier(critical)) < 1e-13
assert abs(critical-3.578234808447358) < 1e-10

print(f"critical_tau={critical:.12f}")
print("PASS: critical threshold located")

for tau in (0.1, 1.0, 3.0, 4.0, 5.0):
    multiplier = exact_multiplier(tau)

    if tau < critical:
        assert multiplier > 1.0
        classification = "UNSTABLE"
    else:
        assert multiplier < 1.0
        classification = "STABLE"

    print(
        f"tau={tau:.2f} "
        f"multiplier={multiplier:.8f} "
        f"classification={classification}"
    )

print("PASS: switching threshold verified")
