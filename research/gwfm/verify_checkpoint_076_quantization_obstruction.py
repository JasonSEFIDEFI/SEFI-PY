#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 076

A polynomial quantization obstruction.

Conditional assumption:
[X,P] = i*hbar*I

Use exact normal-ordered Weyl algebra calculations.
No physical quantum emergence is claimed.
"""

import sympy as sp

h = sp.symbols("hbar", positive=True, real=True)
I = sp.I

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 076")
print("POLYNOMIAL QUANTIZATION OBSTRUCTION")
print("=" * 64)

# Operators are dictionaries:
# (a,b) -> coefficient of X**a P**b.
# All monomials are normal ordered, X before P.

def clean(A):
    return {
        k: sp.simplify(v)
        for k, v in A.items()
        if sp.simplify(v) != 0
    }

def add(A, B, sign=1):
    C = dict(A)
    for k, v in B.items():
        C[k] = C.get(k, 0) + sign*v
    return clean(C)

def scale(A, factor):
    return clean({k: factor*v for k, v in A.items()})

def multiply(A, B):
    C = {}
    for (a,b), va in A.items():
        for (c,d), vb in B.items():
            for k in range(min(b,c)+1):
                coefficient = (
                    va*vb*sp.binomial(b,k)
                    *sp.factorial(c)
                    / sp.factorial(c-k)
                    *(-I*h)**k
                )
                key = (a+c-k,b+d-k)
                C[key] = C.get(key,0)+coefficient
    return clean(C)

def commutator(A,B):
    return add(multiply(A,B),multiply(B,A),-1)

def equal(A,B):
    return clean(add(A,B,-1)) == {}

def monomial(a,b):
    return {(a,b): sp.Integer(1)}

def weyl(a,b):
    # Symmetric/Weyl ordering of x**a * p**b.
    result = {}
    for k in range(min(a,b)+1):
        coefficient = (
            sp.binomial(a,k)
            *sp.binomial(b,k)
            *sp.factorial(k)
            *(-I*h/2)**k
        )
        result[(a-k,b-k)] = coefficient
    return clean(result)

ONE = monomial(0,0)
X = monomial(1,0)
P = monomial(0,1)

# 1. Canonical commutator
assert equal(commutator(X,P),scale(ONE,I*h))
print("PASS: [X,P] = i*hbar*I")

# 2. Classical Poisson identities
x,p = sp.symbols("x p")

def PB(f,g):
    return sp.expand(
        sp.diff(f,x)*sp.diff(g,p)
        -sp.diff(f,p)*sp.diff(g,x)
    )

target = x**2*p**2

assert sp.expand(PB(x**3,p**3)/9-target)==0
assert sp.expand(PB(x**2*p,x*p**2)/3-target)==0

print("PASS: two classical identities for x^2*p^2")

# 3. Weyl ordering is a chosen prescription.
Qtarget = weyl(2,2)

# First bracket route:
# Q_A = [Q(x^3),Q(p^3)]/(9*i*hbar)
QA = scale(
    commutator(weyl(3,0),weyl(0,3)),
    1/(9*I*h)
)

# Second bracket route:
# Q_B = [Q(x^2*p),Q(x*p^2)]/(3*i*hbar)
QB = scale(
    commutator(weyl(2,1),weyl(1,2)),
    1/(3*I*h)
)

print("PASS: both bracket constructions evaluated")

# 4. Verify explicit corrections relative to
# Weyl-ordered x^2*p^2.

assert equal(
    QA,
    add(Qtarget,scale(ONE,-sp.Rational(1,6)*h**2))
)

assert equal(
    QB,
    add(Qtarget,scale(ONE,sp.Rational(1,6)*h**2))
)

print("PASS: first route has -hbar^2/6 correction")
print("PASS: second route has +hbar^2/6 correction")

# 5. The routes disagree.
difference = add(QA,QB,-1)
expected = scale(ONE,-h**2/3)

assert equal(difference,expected)
assert not equal(QA,QB)

print("PASS: routes disagree by -hbar^2/3 * I")

# 6. Classical limit of the discrepancy.
classical_limit = sp.limit(
    difference[(0,0)],h,0
)

assert classical_limit == 0

print("PASS: discrepancy vanishes as hbar -> 0")

# 7. Both routes have the same leading
# classical polynomial term.
assert QA[(2,2)] == 1
assert QB[(2,2)] == 1

print("PASS: both retain the classical leading term")

print("\n" + "=" * 64)
print("CHECKPOINT 076: ALL TESTS PASSED")
print("=" * 64)

print("ESTABLISHED CONDITIONALLY:")
print("- The canonical Weyl algebra is consistent.")
print("- Two classical bracket identities give")
print("  incompatible exact operator assignments.")
print("- The discrepancy is proportional to hbar^2.")
print("- The classical limit removes the discrepancy.")

print("\nNOT ESTABLISHED:")
print("- Quantum mechanics emerges from GWFM.")
print("- The canonical commutator is geometrically forced.")
print("- The value of hbar is derived.")
print("- The Born rule is derived.")
