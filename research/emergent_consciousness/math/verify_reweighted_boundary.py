import itertools
import numpy as np

print("Checkpoint 029: Reweighted metric boundary")

S = np.array([[1.,4.],[0.,1.]])
Sinv = np.linalg.inv(S)

B = S @ np.diag([0.8,0.3]) @ Sinv
C = S @ np.diag([0.4,0.7]) @ Sinv

P = np.array([[1.,-4.],[-4.,41.]])

# T.T @ T = P
T = np.diag([1.,5.]) @ Sinv
Tinv = np.linalg.inv(T)

assert np.allclose(T.T @ T,P)

signs = list(itertools.product([-1.,1.],repeat=4))

def worst_case(eps):
    maximum = -1.
    winner = None

    for name,base in [("B",B),("C",C)]:
        for sign in signs:
            E = eps*np.array(sign).reshape(2,2)
            M = base+E
            sigma = np.linalg.norm(T@M@Tinv,2)

            if sigma > maximum:
                maximum = sigma
                winner = (name,sign)

    return maximum,winner

assert worst_case(0.03)[0] < 1

left,right = 0.03,0.1

while worst_case(right)[0] < 1:
    right *= 2
    assert right < 100

for _ in range(60):
    mid = (left+right)/2

    if worst_case(mid)[0] < 1:
        left = mid
    else:
        right = mid

critical = (left+right)/2
value,winner = worst_case(critical)

assert abs(value-1) < 1e-10

print(f"critical_epsilon={critical:.12f}")
print(f"critical_metric_norm={value:.12f}")
print(f"limiting_operator={winner[0]}")
print(f"limiting_signs={winner[1]}")

for factor in (0.99,1.00,1.01):
    eps = factor*critical
    sigma,_ = worst_case(eps)

    print(
        f"factor={factor:.2f} "
        f"epsilon={eps:.12f} "
        f"worst_sigma={sigma:.12f}"
    )

assert worst_case(0.99*critical)[0] < 1
assert worst_case(1.01*critical)[0] > 1

print("PASS: fixed-metric boundary bracketed")
print("PASS: Checkpoint 029 verified")
