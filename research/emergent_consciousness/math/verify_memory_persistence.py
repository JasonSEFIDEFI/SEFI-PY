from fractions import Fraction
from collections import defaultdict
from math import log2

def mutual_information(joint):
    px = defaultdict(Fraction)
    py = defaultdict(Fraction)

    for (x, y), p in joint.items():
        px[x] += p
        py[y] += p

    result = 0.0
    for (x, y), p in joint.items():
        if p:
            result += float(p) * log2(
                float(p / (px[x] * py[y]))
            )
    return result

def measure(step):
    joint = defaultdict(Fraction)

    for e0 in (0, 1):
        m = 0
        e = e0

        for _ in range(step):
            m = m ^ e

        joint[(e0, m)] += Fraction(1, 2)

    return mutual_information(joint)

for step in range(5):
    print(
        f"Step {step}: "
        f"I(E0;M{step}) = {measure(step):.6f} bits"
    )

assert abs(measure(1) - 1.0) < 1e-12
assert abs(measure(2)) < 1e-12
assert abs(measure(3) - 1.0) < 1e-12

print("PASS: reversible memory oscillation verified")
