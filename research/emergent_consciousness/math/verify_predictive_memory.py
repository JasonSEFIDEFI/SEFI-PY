from fractions import Fraction
from math import log2
from collections import defaultdict

joint = defaultdict(Fraction)

for u0 in (0, 1):
    for v0 in (0, 1):
        u1 = v0
        u2 = u0
        joint[(u0, u1, u2)] += Fraction(1, 4)

def entropy(indices):
    probabilities = defaultdict(Fraction)
    for state, p in joint.items():
        key = tuple(state[i] for i in indices)
        probabilities[key] += p
    return -sum(float(p) * log2(float(p))
                for p in probabilities.values())

# I(U2;U1)
mi_current = (
    entropy((2,)) + entropy((1,))
    - entropy((1, 2))
)

# I(U2;U0 | U1)
cmi_past = (
    entropy((2, 1)) + entropy((0, 1))
    - entropy((1,)) - entropy((0, 1, 2))
)

print("Predictive memory verification")
print(f"I(U2; U1) = {mi_current:.6f} bits")
print(f"I(U2; U0 | U1) = {cmi_past:.6f} bits")

assert abs(mi_current) < 1e-12
assert abs(cmi_past - 1.0) < 1e-12

print("PASS: predictive information verified")
