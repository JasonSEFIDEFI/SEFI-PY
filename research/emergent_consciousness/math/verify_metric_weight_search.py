import itertools
import numpy as np

print("Checkpoint 032: Metric weight search")

S = np.array([[1.,4.],[0.,1.]])
Sinv = np.linalg.inv(S)

B = S @ np.diag([0.8,0.3]) @ Sinv
C = S @ np.diag([0.4,0.7]) @ Sinv

signs = np.array(
    list(itertools.product([-1.,1.],repeat=4))
).reshape(-1,2,2)

def worst_norm(eps,w):
    T = np.diag([1.,np.sqrt(w)]) @ Sinv
    Tinv = np.linalg.inv(T)

    worst = 0.

    for base in (B,C):
        for pattern in signs:
            M = base + eps*pattern
            value = np.linalg.norm(T@M@Tinv,2)
            worst = max(worst,value)

    return worst

# Reproduce previously verified cases.
assert worst_norm(0.03,25.) < 1
assert worst_norm(0.035,25.) > 1

print("PASS: previous metric results recovered")

# Search a broad, finite family of weights.
weights = np.geomspace(0.1,10000.,301)

def best_at(eps):
    values = [
        worst_norm(eps,w) for w in weights
    ]
    index = int(np.argmin(values))
    return weights[index],values[index]

for eps in (0.03,0.035,0.04,0.045,0.05):
    w,value = best_at(eps)
    print(
        f"epsilon={eps:.6f} "
        f"best_grid_weight={w:.9f} "
        f"worst_norm={value:.12f} "
        f"certified={value < 1}"
    )

print("PASS: restricted grid search completed")
print("Checkpoint 032: exploratory optimization complete")
