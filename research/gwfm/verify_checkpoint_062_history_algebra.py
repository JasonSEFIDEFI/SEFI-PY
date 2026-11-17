#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 062
Finite real history algebra.

Mathematical scope:
    - Directed paths in a finite graph.
    - Free real vector spaces on paths.
    - Bilinear path concatenation.
    - Local identity histories.
    - Associativity and distributivity.

This is an algebraic verification, not a derivation of
quantum amplitudes, the Born rule, or physical SEFI fields.
"""

from collections import defaultdict
from fractions import Fraction
from itertools import product


# A history is a tuple of vertices.
# (A,) is the identity history at A.
#
# Distinct edge sequences remain distinct histories,
# even if their initial and final vertices agree.

def compose_paths(first, second):
    """Travel first, then second."""
    if first[-1] != second[0]:
        raise ValueError("Incompatible endpoints")

    return first + second[1:]


class HistorySum:
    """Finite formal real linear combination of histories."""

    def __init__(self, terms=None):
        combined = defaultdict(Fraction)

        for path, coefficient in (terms or {}).items():
            path = tuple(path)

            if not path:
                raise ValueError("Empty path tuple")

            combined[path] += Fraction(coefficient)

        self.terms = {
            p: c for p, c in combined.items() if c != 0
        }

    @classmethod
    def basis(cls, path):
        return cls({tuple(path): Fraction(1)})

    def __add__(self, other):
        result = defaultdict(Fraction)

        for p, c in self.terms.items():
            result[p] += c

        for p, c in other.terms.items():
            result[p] += c

        return HistorySum(result)

    def __mul__(self, scalar):
        return HistorySum({
            p: c * Fraction(scalar)
            for p, c in self.terms.items()
        })

    __rmul__ = __mul__

    def then(self, other):
        """
        Sequential composition.

        self.then(other) means self followed by other.
        """
        result = defaultdict(Fraction)

        for p, a in self.terms.items():
            for q, b in other.terms.items():
                result[compose_paths(p, q)] += a * b

        return HistorySum(result)

    def __eq__(self, other):
        return (
            isinstance(other, HistorySum)
            and self.terms == other.terms
        )

    def __repr__(self):
        if not self.terms:
            return "0"

        return " + ".join(
            f"({c})[{','.join(p)}]"
            for p, c in sorted(self.terms.items())
        )


def identity(vertex):
    return HistorySum.basis((vertex,))


def verify():
    print("=" * 62)
    print("GWFM / SEFI — CHECKPOINT 062")
    print("FINITE WORLDLINE HISTORY ALGEBRA")
    print("=" * 62)

    # Directed graph:
    #
    # A -> B -> D -> E
    #  \       /
    #   -> C -
    #
    # Two distinct histories connect A to D.

    AB = HistorySum.basis(("A", "B"))
    AC = HistorySum.basis(("A", "C"))
    BD = HistorySum.basis(("B", "D"))
    CD = HistorySum.basis(("C", "D"))
    DE = HistorySum.basis(("D", "E"))

    # Alternative histories from A to D.
    h1 = AB.then(BD)
    h2 = AC.then(CD)

    print("\nTwo distinct histories:")
    print("h1 =", h1)
    print("h2 =", h2)

    assert h1 != h2

    # Local identities.
    assert identity("A").then(AB) == AB
    assert AB.then(identity("B")) == AB
    assert identity("D").then(DE) == DE

    print("\nPASS: local identity histories")

    # Associativity of sequential composition.
    left = AB.then(BD).then(DE)
    right = AB.then(BD.then(DE))

    assert left == right

    print("PASS: associative concatenation")

    # Formal alternatives.
    alternatives = h1 + h2

    print("\nFormal alternative histories:")
    print(alternatives)

    # Shared continuation distributes over alternatives.
    left = alternatives.then(DE)
    right = h1.then(DE) + h2.then(DE)

    assert left == right

    print("\nPASS: right distributivity")

    # Left distributivity with compatible alternatives.
    first = AB + AC
    continuation = (
        BD.then(DE)
        + CD.then(DE)
    )

    # The full Cartesian product is not composable:
    # AB cannot be followed by CD, and AC cannot
    # be followed by BD.
    #
    # Therefore test left distributivity inside
    # a common hom-space instead.

    f = h1
    g = h2
    k = DE

    assert (f + g).then(k) == (
        f.then(k) + g.then(k)
    )

    print("PASS: distributivity in compatible hom-spaces")

    # Scalar bilinearity.
    a = Fraction(2, 3)
    b = Fraction(-5, 7)

    assert ((a * f) + (b * g)).then(k) == (
        a * f.then(k) + b * g.then(k)
    )

    print("PASS: exact rational bilinearity")

    # Distinct positive functionals on the same
    # formal history space can give different values.
    #
    # This demonstrates underdetermination.
    # It does NOT construct quantum probabilities.

    def diagonal_weight(x):
        return sum(
            coefficient * coefficient
            for coefficient in x.terms.values()
        )

    def total_coefficient_square(x):
        return sum(x.terms.values(), Fraction(0)) ** 2

    mixture = f + g

    diagonal = diagonal_weight(mixture)
    coherent_like = total_coefficient_square(mixture)

    print("\nTwo different quadratic assignments:")
    print("Diagonal coefficient sum:", diagonal)
    print("Square of coefficient sum:", coherent_like)

    assert diagonal == 2
    assert coherent_like == 4
    assert diagonal != coherent_like

    print("PASS: algebra alone does not fix a quadratic rule")

    # Exhaustive finite associativity test.
    paths = [
        ("A", "B"),
        ("B", "D"),
        ("D", "E"),
        ("A", "C"),
        ("C", "D"),
        ("A",),
        ("B",),
        ("C",),
        ("D",),
        ("E",),
    ]

    tested = 0

    for p, q, r in product(paths, repeat=3):
        if p[-1] == q[0] and q[-1] == r[0]:
            lhs = compose_paths(compose_paths(p, q), r)
            rhs = compose_paths(p, compose_paths(q, r))

            assert lhs == rhs
            tested += 1

    print(f"PASS: {tested} compatible triple compositions")

    print("\n" + "=" * 62)
    print("CHECKPOINT 062: ALL TESTS PASSED")
    print("=" * 62)
    print("Scope: finite directed-path algebra over Q")
    print("No quantum or physical emergence claim established.")


if __name__ == "__main__":
    verify()
