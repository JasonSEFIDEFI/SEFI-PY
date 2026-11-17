#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 065
Positive linear functionals on a formal path *-algebra.

Scope:
- Finite directed graph with formally adjoined reverse edges
- Free rational path algebra
- Two positive normalized linear functionals
- Exact arithmetic and non-uniqueness tests

No physical emergence or Born-rule derivation is claimed.
"""

from collections import defaultdict
from fractions import Fraction as Q
from itertools import product

V = ("A", "B", "C", "D")

FORWARD = {
    ("A", "B"), ("A", "C"),
    ("B", "D"), ("C", "D")
}

EDGES = FORWARD | {(v, u) for u, v in FORWARD}


def valid(p):
    return bool(p) and all(
        (p[i], p[i + 1]) in EDGES
        for i in range(len(p) - 1)
    )


class Algebra:
    def __init__(self, terms=None):
        d = defaultdict(Q)
        for p, c in (terms or {}).items():
            p = tuple(p)
            if not valid(p):
                raise ValueError(f"Invalid path: {p}")
            d[p] += Q(c)
        self.terms = {p: c for p, c in d.items() if c}

    @staticmethod
    def basis(path):
        return Algebra({tuple(path): 1})

    def __add__(self, other):
        d = defaultdict(Q)
        for a in (self, other):
            for p, c in a.terms.items():
                d[p] += c
        return Algebra(d)

    def __rmul__(self, scalar):
        return Algebra({
            p: Q(scalar) * c
            for p, c in self.terms.items()
        })

    def __sub__(self, other):
        return self + (-1) * other

    def then(self, other):
        d = defaultdict(Q)
        for p, a in self.terms.items():
            for q, b in other.terms.items():
                if p[-1] == q[0]:
                    d[p + q[1:]] += a * b
        return Algebra(d)

    def star(self):
        return Algebra({
            p[::-1]: c for p, c in self.terms.items()
        })

    def __eq__(self, other):
        return (
            isinstance(other, Algebra)
            and self.terms == other.terms
        )


one = Algebra({(v,): 1 for v in V})


def omega_vertex(a, vertex="A"):
    """Coefficient of the length-zero identity at vertex."""
    return a.terms.get((vertex,), Q(0))


def endpoint_matrix(a):
    """
    Represent each path p as E_(start,end).

    Distinct paths with identical endpoints have
    the same image under this representation.
    """
    M = {(u, v): Q(0) for u in V for v in V}
    for p, c in a.terms.items():
        M[(p[0], p[-1])] += c
    return M


WEIGHTS = {v: Q(1, 2) for v in V}


def omega_endpoint(a):
    """Vector state <w, pi(a) w>, with ||w||^2 = 1."""
    M = endpoint_matrix(a)
    return sum(
        (
            WEIGHTS[u] * M[(u, v)] * WEIGHTS[v]
            for u in V for v in V
        ),
        Q(0)
    )


def quadratic(a, omega):
    return omega(a.star().then(a))


def main():
    print("=" * 64)
    print("GWFM / SEFI — CHECKPOINT 065")
    print("POSITIVE LINEAR FUNCTIONALS")
    print("=" * 64)

    assert omega_vertex(one) == 1
    assert omega_endpoint(one) == 1
    print("PASS: both functionals normalized")

    h1 = Algebra.basis(("A", "B", "D"))
    h2 = Algebra.basis(("A", "C", "D"))

    examples = [
        Algebra(),
        one,
        h1,
        h2,
        h1 + h2,
        h1 - h2,
        2 * h1 - 3 * h2,
        Algebra.basis(("B", "A")) + h1,
    ]

    for a in examples:
        assert quadratic(a, omega_vertex) >= 0
        assert quadratic(a, omega_endpoint) >= 0

    print("PASS: positivity on selected algebra elements")

    # Check linearity exactly on sample elements.
    for omega in (omega_vertex, omega_endpoint):
        for a, b in product(examples, repeat=2):
            assert omega(a + b) == omega(a) + omega(b)
            assert omega(3 * a) == 3 * omega(a)

    print("PASS: exact rational linearity tests")

    # The endpoint representation respects
    # multiplication and the involution.
    for a, b in product(examples, repeat=2):
        A = endpoint_matrix(a)
        B = endpoint_matrix(b)
        AB = endpoint_matrix(a.then(b))

        for u, v in product(V, repeat=2):
            expected = sum(
                (A[(u, k)] * B[(k, v)] for k in V),
                Q(0)
            )
            assert AB[(u, v)] == expected

        AS = endpoint_matrix(a.star())
        for u, v in product(V, repeat=2):
            assert AS[(u, v)] == A[(v, u)]

    print("PASS: endpoint map respects product and star")

    # Both histories have the same endpoints.
    assert h1 != h2
    assert endpoint_matrix(h1) == endpoint_matrix(h2)

    print("PASS: distinct paths share endpoint image")

    print("\nQuadratic values omega(a* a):")
    print("                       Vertex     Endpoint")

    for label, a in [
        ("h1", h1),
        ("h2", h2),
        ("h1 + h2", h1 + h2),
        ("h1 - h2", h1 - h2),
    ]:
        v = quadratic(a, omega_vertex)
        e = quadratic(a, omega_endpoint)
        print(f"{label:<22} {str(v):<10} {e}")

    assert quadratic(h1, omega_vertex) == 0
    assert quadratic(h1, omega_endpoint) == Q(1, 4)
    assert quadratic(h1 + h2, omega_endpoint) == 1
    assert quadratic(h1 - h2, omega_endpoint) == 0

    print("\nPASS: positive states are non-unique")
    print("PASS: endpoint state has coherent cross terms")

    # An exact algebraic formula proves positivity
    # for arbitrary finite rational path sums.
    #
    # For the vertex state:
    # omega_A(a* a) = coefficient(a, e_A)^2.
    #
    # For the endpoint vector state:
    # omega_w(a* a) = sum_v |(pi(a)w)_v|^2.
    #
    # These identities are valid beyond tested samples.

    print("\nTHEOREM (conditional on formal extension):")
    print("Both constructions are positive normalized")
    print("linear functionals on the rational path *-algebra.")

    print("\nLIMITATIONS:")
    print("- Reverse edges were introduced formally.")
    print("- Endpoint representation was selected.")
    print("- Neither state is derived from spacetime geometry.")
    print("- No physical interference or Born rule is derived.")

    print("\n" + "=" * 64)
    print("CHECKPOINT 065: ALL TESTS PASSED")
    print("=" * 64)


if __name__ == "__main__":
    main()
