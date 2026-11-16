from fractions import Fraction
from math import log2

def entropy_binary(q):
    if q == 0 or q == 1:
        return 0.0
    x = float(q)
    return -x*log2(x) - (1-x)*log2(1-x)

def iterative_error(p, n):
    q = Fraction(0)
    for _ in range(n):
        q = p + (1 - 2*p)*q
    return q

def analytical_error(p, n):
    return (1 - (1 - 2*p)**n) / 2

print("Checkpoint 009: Repeated independent noise")

for p in (
    Fraction(0),
    Fraction(1, 10),
    Fraction(1, 2),
    Fraction(9, 10),
    Fraction(1)
):
    for n in (0, 1, 2, 5, 10):
        q_iter = iterative_error(p, n)
        q_exact = analytical_error(p, n)

        assert q_iter == q_exact

        info = 1 - entropy_binary(q_exact)

        if p in (Fraction(1, 10), Fraction(1)):
            print(
                f"p={float(p):.1f}, n={n:2d}, "
                f"q={float(q_exact):.6f}, "
                f"I={info:.6f} bits"
            )

print("PASS: repeated-noise recurrence verified")
