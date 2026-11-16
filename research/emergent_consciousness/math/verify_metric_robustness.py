import numpy as np
import sympy as sp

print("Checkpoint 025: Metric robustness")

S = sp.Matrix([[1, 4], [0, 1]])
P = S.inv().T * S.inv()

D1 = sp.diag(sp.Rational(4, 5), sp.Rational(3, 10))
D2 = sp.diag(sp.Rational(2, 5), sp.Rational(7, 10))

B = S * D1 * S.inv()
C = S * D2 * S.inv()

eps = sp.Rational(1, 200)
kappa = 9 + 4 * sp.sqrt(5)
r = sp.Rational(4, 5) + 2 * kappa * eps
q = sp.simplify(r**2)

assert r < 1
assert P.det() > 0

print("PASS: analytical robustness bound")
print(f"condition_number={float(kappa):.9f}")
print(f"norm_bound={float(r):.9f}")
print(f"energy_bound={float(q):.9f}")

# Independent numerical perturbation tests
rng = np.random.default_rng(20261115)

Sn = np.array(S, dtype=float)
Sinv = np.linalg.inv(Sn)

bases = [
    np.array(B, dtype=float),
    np.array(C, dtype=float)
]

maximum = 0.0

for _ in range(10000):
    base = bases[int(rng.integers(0, 2))]
    E = rng.uniform(-float(eps), float(eps), (2, 2))

    transformed = Sinv @ (base + E) @ Sn
    ratio = np.linalg.norm(transformed, 2)**2

    maximum = max(maximum, ratio)
    assert ratio <= float(q) + 1e-12

print(f"sampled_max_energy_ratio={maximum:.9f}")
print("PASS: independent numerical tests")
print("PASS: Checkpoint 025 verified")
