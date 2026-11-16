from fractions import Fraction
from collections import defaultdict
from math import log2

def mutual_information(joint):
    px = defaultdict(Fraction)
    py = defaultdict(Fraction)

    for (x, y), p in joint.items():
        px[x] += p
        py[y] += p

    total = 0.0
    for (x, y), p in joint.items():
        if p:
            total += float(p) * log2(
                float(p / (px[x] * py[y]))
            )
    return total

def experiment(memory_values):
    joint = defaultdict(Fraction)
    p = Fraction(1, 2 * len(memory_values))

    for e0 in (0, 1):
        for m0 in memory_values:
            m1 = m0 ^ e0
            joint[(e0, m1)] += p

    return mutual_information(joint)

random_memory = experiment((0, 1))
initialized_memory = experiment((0,))

print("Information acquisition verification")
print(f"Random initial memory: {random_memory:.6f} bits")
print(f"Initialized memory: {initialized_memory:.6f} bits")

assert abs(random_memory) < 1e-12
assert abs(initialized_memory - 1.0) < 1e-12

print("PASS: initial-condition dependence verified")
