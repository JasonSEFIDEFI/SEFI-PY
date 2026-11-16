from itertools import product
from fractions import Fraction

def majority(bits):
    return int(sum(bits) >= 2)

def exact_error(p):
    error = Fraction(0)

    for noise in product((0, 1), repeat=3):
        flips = sum(noise)

        probability = (
            p**flips * (1-p)**(3-flips)
        )

        original = 0
        received = tuple(
            original ^ bit for bit in noise
        )

        if majority(received) != original:
            error += probability

    return error

print("Checkpoint 010: Redundant memory")

for p in (
    Fraction(0),
    Fraction(1, 10),
    Fraction(1, 4),
    Fraction(1, 2),
    Fraction(3, 4),
    Fraction(1)
):
    observed = exact_error(p)
    predicted = 3*p*p - 2*p*p*p

    assert observed == predicted

    print(
        f"p={float(p):.2f}, "
        f"error={float(observed):.6f}"
    )

print("PASS: repetition-code identity verified")
