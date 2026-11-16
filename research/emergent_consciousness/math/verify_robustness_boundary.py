import itertools
import numpy as np

print("Checkpoint 026: Exact box-vertex robustness")

S = np.array([[1., 4.], [0., 1.]])
Sinv = np.linalg.inv(S)

D1 = np.diag([0.8, 0.3])
D2 = np.diag([0.4, 0.7])

B = S @ D1 @ Sinv
C = S @ D2 @ Sinv

signs = list(itertools.product([-1., 1.], repeat=4))

def worst_case(eps):
    maximum = -1.
    winner = None

    for name, base in [("B", B), ("C", C)]:
        for sign in signs:
            E = eps * np.array(sign).reshape(2, 2)
            A = Sinv @ (base + E) @ S
            sigma = np.linalg.norm(A, 2)

            if sigma > maximum:
                maximum = sigma
                winner = (name, sign)

    return maximum, winner

kappa = np.linalg.cond(S, 2)
conservative = 0.2 / (2 * kappa)

assert np.isclose(
    kappa, 9 + 4*np.sqrt(5), rtol=1e-12
)

assert np.isclose(worst_case(0)[0], 0.8)
assert worst_case(conservative)[0] <= 1 + 1e-12

print(f"conservative_epsilon={conservative:.12f}")
print("PASS: conservative threshold verified")

# Bracket the first loss of strict contraction.
left = 0.
right = 0.1

assert worst_case(left)[0] < 1

while worst_case(right)[0] <= 1:
    right *= 2
    assert right < 100

for _ in range(60):
    mid = (left + right)/2

    if worst_case(mid)[0] < 1:
        left = mid
    else:
        right = mid

critical = (left + right)/2
sigma, winner = worst_case(critical)

assert abs(sigma - 1) < 1e-10
assert critical >= conservative

print(f"critical_epsilon={critical:.12f}")
print(f"critical_sigma={sigma:.12f}")
print(f"limiting_operator={winner[0]}")
print(f"limiting_signs={winner[1]}")
print("PASS: worst-case boundary located")

for factor in [0.99, 1.00, 1.01]:
    eps = factor * critical
    value, _ = worst_case(eps)

    print(
        f"factor={factor:.2f} "
        f"epsilon={eps:.12f} "
        f"worst_sigma={value:.12f}"
    )

assert worst_case(0.99*critical)[0] < 1
assert worst_case(1.01*critical)[0] > 1

print("PASS: transition verified")

# Explicitly verify the finite vertex count.
assert len(signs) * 2 == 32

print("PASS: all 32 vertex families considered")
print("PASS: Checkpoint 026 verified")
