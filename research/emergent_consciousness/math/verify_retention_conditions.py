from fractions import Fraction
from collections import defaultdict
from math import log2

def entropy(distribution):
    return -sum(
        float(p) * log2(float(p))
        for p in distribution.values() if p
    )

def information_for(rule):
    joint = defaultdict(Fraction)

    for e0 in (0, 1):
        e1 = e0
        m1 = e0

        if rule == "identity":
            m2 = m1
        elif rule == "repeat_xor":
            m2 = m1 ^ e1
        elif rule == "flip":
            m2 = m1 ^ 1
        else:
            raise ValueError(rule)

        joint[(e0, m2)] += Fraction(1, 2)

    pe = defaultdict(Fraction)
    pm = defaultdict(Fraction)

    for (e, m), p in joint.items():
        pe[e] += p
        pm[m] += p

    h_e = entropy(pe)
    h_m = entropy(pm)
    h_em = entropy(joint)

    return h_e + h_m - h_em, h_em - h_m

expected = {
    "identity": (1.0, 0.0),
    "repeat_xor": (0.0, 1.0),
    "flip": (1.0, 0.0),
}

print("Checkpoint 006: Retention conditions")

for rule, (expected_i, expected_h) in expected.items():
    info, conditional = information_for(rule)

    print(
        f"{rule}: "
        f"I(E0;M2)={info:.6f} bits, "
        f"H(E0|M2)={conditional:.6f} bits"
    )

    assert abs(info - expected_i) < 1e-12
    assert abs(conditional - expected_h) < 1e-12

print("PASS: information retention conditions verified")
