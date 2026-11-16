from itertools import product
from fractions import Fraction
from collections import defaultdict
from math import log2

def entropy(probs):
    return -sum(
        float(p) * log2(float(p))
        for p in probs if p
    )

def conditional_entropy(encoding):
    joint = defaultdict(Fraction)
    memory = defaultdict(Fraction)

    for e in (0, 1):
        m = encoding[e]
        joint[(e, m)] += Fraction(1, 2)
        memory[m] += Fraction(1, 2)

    return entropy(joint.values()) - entropy(memory.values())

def has_decoder(encoding):
    # Both environmental values must be recoverable.
    return encoding[0] != encoding[1]

count = 0
perfect = 0

for encoding in product(range(3), repeat=2):
    h = conditional_entropy(encoding)
    decoder_exists = has_decoder(encoding)

    assert (abs(h) < 1e-12) == decoder_exists

    count += 1
    perfect += int(decoder_exists)

print("Checkpoint 007: General retention theorem")
print(f"Encodings examined: {count}")
print(f"Perfectly recoverable: {perfect}")
print(f"Not recoverable: {count-perfect}")
print("PASS: decoder equivalence verified")
