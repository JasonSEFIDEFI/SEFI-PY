#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 064

Formal path reversal and involution.

Scope:
Finite directed graph, free path algebra over Q,
and its formal reverse-edge extension.

No physical time reversal, antiparticle interpretation,
or quantum emergence is assumed.
"""

from collections import defaultdict
from fractions import Fraction


VERTICES = ("A", "B", "C", "D")

FORWARD_EDGES = {
    ("A", "B"),
    ("A", "C"),
    ("B", "D"),
    ("C", "D"),
}

REVERSE_EDGES = {
    (v, u) for u, v in FORWARD_EDGES
}

EXTENDED_EDGES = FORWARD_EDGES | REVERSE_EDGES


def valid_path(path, edges):
    if not path:
        return False

    return all(
        (path[i], path[i + 1]) in edges
        for i in range(len(path) - 1)
    )


def reverse_path(path):
    return tuple(reversed(path))


def compose_paths(p, q):
    if p[-1] != q[0]:
        raise ValueError("Incompatible endpoints")

    return p + q[1:]


class PathAlgebra:
    def __init__(self, terms=None):
        coefficients = defaultdict(Fraction)

        for path, coefficient in (terms or {}).items():
            path = tuple(path)

            if not valid_path(path, EXTENDED_EDGES):
                raise ValueError(f"Invalid path: {path}")

            coefficients[path] += Fraction(coefficient)

        self.terms = {
            p: c for p, c in coefficients.items()
            if c != 0
        }

    @classmethod
    def basis(cls, path):
        return cls({tuple(path): 1})

    def __add__(self, other):
        result = defaultdict(Fraction)

        for p, c in self.terms.items():
            result[p] += c

        for p, c in other.terms.items():
            result[p] += c

        return PathAlgebra(result)

    def __mul__(self, scalar):
        return PathAlgebra({
            p: c * Fraction(scalar)
            for p, c in self.terms.items()
        })

    __rmul__ = __mul__

    def then(self, other):
        result = defaultdict(Fraction)

        for p, a in self.terms.items():
            for q, b in other.terms.items():
                if p[-1] != q[0]:
                    continue

                result[compose_paths(p, q)] += a * b

        return PathAlgebra(result)

    def star(self):
        return PathAlgebra({
            reverse_path(p): c
            for p, c in self.terms.items()
        })

    def __eq__(self, other):
        return (
            isinstance(other, PathAlgebra)
            and self.terms == other.terms
        )


def verify():
    print("=" * 64)
    print("GWFM / SEFI — CHECKPOINT 064")
    print("FORMAL HISTORY INVOLUTION")
    print("=" * 64)

    p = ("A", "B", "D")

    assert valid_path(p, FORWARD_EDGES)
    assert not valid_path(reverse_path(p), FORWARD_EDGES)

    print("PASS: original forward paths not closed under reversal")

    assert valid_path(reverse_path(p), EXTENDED_EDGES)

    print("PASS: formal reverse-edge extension closes reversal")

    a = PathAlgebra.basis(("A", "B"))
    b = PathAlgebra.basis(("B", "D"))

    assert a.star().star() == a

    print("PASS: involution squared is identity")

    left = a.then(b).star()
    right = b.star().then(a.star())

    assert left == right

    print("PASS: involution reverses composition order")

    x = a + 2 * PathAlgebra.basis(("A", "C"))
    y = b + 3 * PathAlgebra.basis(("C", "D"))

    assert (x + y).star() == x.star() + y.star()

    print("PASS: involution is additive")

    assert (x.then(y)).star() == y.star().then(x.star())

    print("PASS: involution respects formal path sums")

    # A simple positive functional is obtained from
    # the left regular representation on finite paths.
    #
    # We do not construct that infinite-dimensional
    # representation here. Instead, verify a finite
    # coefficient norm as an illustrative nonnegative
    # quadratic expression.

    def coefficient_norm_squared(element):
        return sum(
            c * c for c in element.terms.values()
        )

    assert coefficient_norm_squared(x) >= 0

    print("PASS: coefficient norm is nonnegative")
    print("NOTE: coefficient norm alone is not a")
    print("      constructed positive linear functional.")

    print("\nCHECKPOINT 064: ALL TESTS PASSED")
    print("SCOPE: formal reverse-edge path algebra over Q")
    print("No physical interpretation established.")


if __name__ == "__main__":
    verify()
