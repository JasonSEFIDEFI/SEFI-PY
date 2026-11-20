import sympy as S
from time import perf_counter

start = perf_counter()

print("=== CP098-T90C-R28B ===", flush=True)

W = S.symbols("W", positive=True)

numerator = (
    476656*W**3
    + 20196757*W**2
    + 41330688*W
    + 531948096
)

denominator = 14776336*W

D = S.cancel(numerator/denominator)

print("STAGE 1: exact R28A regression")

expected = {
    S.Rational(1, 4):
        S.Rational(8696808209, 59105344),
    S.Integer(1):
        S.Rational(593952197, 14776336),
    S.Integer(4):
        S.Rational(65682809, 3694084)
}

for value, target in expected.items():
    actual = S.cancel(D.subs(W, value))
    assert actual == target, (value, actual, target)
    print("PASS: W =", value)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: positivity certificate")

poly = S.Poly(numerator, W)

assert all(c > 0 for c in poly.all_coeffs())
assert S.ask(S.Q.positive(W)) is True
assert S.ask(S.Q.positive(denominator)) is True

terms = [
    476656*W**3,
    20196757*W**2,
    41330688*W,
    S.Integer(531948096)
]

assert S.expand(sum(terms) - numerator) == 0

for term in terms:
    assert S.ask(S.Q.positive(term)) is True

print("PASS: numerator positive for W>0")
print("PASS: denominator positive for W>0")
print("PASS: D(W)>0 for all W>0")
print("STAGE 2 PASS", flush=True)

print("STAGE 3: strict lower bound")

lower_bound = S.Rational(
    41330688, 14776336
)

difference = S.cancel(D - lower_bound)

expected_difference = (
    476656*W**2
    + 20196757*W
    + 531948096/W
)/14776336

assert S.cancel(
    difference - expected_difference
) == 0

assert S.ask(
    S.Q.positive(expected_difference)
) is True

print("D(W) =", D)
print("lower bound =", lower_bound)
print("PASS: D(W)>lower bound>0")
print("STAGE 3 PASS", flush=True)

print()
print("R28B COMPLETE")
print("R28A EXPRESSION: REGRESSION VERIFIED")
print("GLOBAL W-FAMILY POSITIVITY: PROVED")
print("D=0 IN THIS FAMILY: IMPOSSIBLE")
print("JOINT W,P² ANALYSIS: PENDING")
print("GENERAL IA BRANCH: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter() - start, 3)
)
