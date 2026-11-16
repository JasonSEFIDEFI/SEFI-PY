from fractions import Fraction
from math import log2

def entropy(probabilities):
    return -sum(
        float(q) * log2(float(q))
        for q in probabilities if q > 0
    )

def exact_information(p):
    # E is uniform; N is independent Bernoulli(p).
    joint = {}
    for e in (0, 1):
        for n, pn in ((0, 1-p), (1, p)):
            m = e ^ n
            joint[(e, m)] = Fraction(1, 2) * pn

    marginal_e = [
        sum(q for (e, m), q in joint.items() if e == x)
        for x in (0, 1)
    ]
    marginal_m = [
        sum(q for (e, m), q in joint.items() if m == x)
        for x in (0, 1)
    ]

    return (
        entropy(marginal_e)
        + entropy(marginal_m)
        - entropy(joint.values())
    )

def binary_entropy(p):
    return entropy((p, 1-p))

print("Checkpoint 008: Noisy information retention")

for p in (
    Fraction(0), Fraction(1, 10),
    Fraction(1, 4), Fraction(1, 2),
    Fraction(3, 4), Fraction(9, 10),
    Fraction(1)
):
    measured = exact_information(p)
    predicted = 1 - binary_entropy(p)

    assert abs(measured - predicted) < 1e-12

    print(
        f"p={float(p):.2f} "
        f"I={measured:.6f} bits"
    )

print("PASS: noisy retention identity verified")
