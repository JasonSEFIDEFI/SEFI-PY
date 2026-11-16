import numpy as np

print("Checkpoint 027: Certificate versus stability")

S = np.array([[1., 4.], [0., 1.]])
Sinv = np.linalg.inv(S)
P = Sinv.T @ Sinv

D = np.diag([0.8, 0.3])
B = S @ D @ Sinv

critical = 0.017967017631
eps = 1.01 * critical

E = eps * np.array([
    [1., 1.],
    [-1., -1.]
])

M = B + E

# Check failure of the original metric.
old_difference = M.T @ P @ M - P
old_max = np.linalg.eigvalsh(old_difference)[-1]

assert old_max > 0

print(f"old_metric_max_eigenvalue={old_max:.12f}")
print("PASS: original strict certificate fails")

# Independent spectral stability test.
eigenvalues = np.linalg.eigvals(M)
rho = max(abs(eigenvalues))

print(f"spectral_radius={rho:.12f}")

assert rho < 1

print("PASS: perturbed operator remains stable")

# Solve M.T Q M - Q = -I.
# vec(M.T Q M) = kron(M.T,M.T) vec(Q)
# using column-major vectorization.
I = np.eye(2)

A = np.eye(4) - np.kron(M.T, M.T)
rhs = I.reshape(4, order="F")

Q = np.linalg.solve(A, rhs)
Q = Q.reshape((2, 2), order="F")
Q = (Q + Q.T) / 2

assert np.linalg.eigvalsh(Q)[0] > 0

residual = M.T @ Q @ M - Q + I
error = np.linalg.norm(residual, ord="fro")

assert error < 1e-10

print(f"new_metric_min_eigenvalue="
      f"{np.linalg.eigvalsh(Q)[0]:.12f}")
print(f"lyapunov_residual={error:.3e}")

print("PASS: new positive-definite metric")
print("PASS: discrete Lyapunov identity")

# Independent finite-series cross-check.
Qseries = np.zeros((2, 2))
power = np.eye(2)

for _ in range(1000):
    Qseries += power.T @ power
    power = M @ power

assert np.linalg.norm(
    Qseries - Q, ord="fro"
) < 1e-8

print("PASS: independent Lyapunov series")
print("PASS: Checkpoint 027 verified")
