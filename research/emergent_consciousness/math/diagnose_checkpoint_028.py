import itertools
import numpy as np

print("Checkpoint 028: Candidate metric diagnosis")

S = np.array([[1., 4.], [0., 1.]])
Sinv = np.linalg.inv(S)

B = S @ np.diag([0.8, 0.3]) @ Sinv
C = S @ np.diag([0.4, 0.7]) @ Sinv

P = np.array([[1., -4.], [-4., 41.]])
eps = 0.03

assert np.min(np.linalg.eigvalsh(P)) > 0

worst = (float("inf"), None, None)
passed = 0

for name, base in [("B", B), ("C", C)]:
    for signs in itertools.product([-1., 1.], repeat=4):
        E = eps * np.array(signs).reshape(2, 2)
        M = base + E

        H = P - M.T @ P @ M
        smallest = np.min(np.linalg.eigvalsh(H))

        if smallest > 0:
            passed += 1

        if smallest < worst[0]:
            worst = (smallest, name, signs)

print(f"positive_definite_vertices={passed}/32")
print(f"worst_min_eigenvalue={worst[0]:.12f}")
print(f"worst_operator={worst[1]}")
print(f"worst_sign_pattern={worst[2]}")

if passed == 32:
    print("PASS: Candidate common metric verified")
else:
    print("NOT VERIFIED: Candidate metric fails")
    print("Checkpoint 028 remains open")
