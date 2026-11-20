import sympy as S

print("=== CP098-T90C-R29I ===")

A = S.Rational(
    2196291909748885,
    111262356920464396
)
B = S.Rational(
    1152303896300159,
    9250279092157
)

X = S.Integer(144209559530)
Y = S.Integer(194)

assert A > 0 and B > 0
assert X > 0 and Y > 0

difference = A**2 * X - B**2 * Y

assert difference > 0
print("PASS: exact squared inequality")

plus = -A*S.sqrt(X) + B*S.sqrt(Y)
minus = -A*S.sqrt(X) - B*S.sqrt(Y)

assert S.simplify(plus) != 0
assert S.simplify(minus) != 0

assert S.N(plus, 30) < 0
assert S.N(minus, 30) < 0

print("PASS: both branches strictly negative")
print("PLUS =", S.N(plus, 20))
print("MINUS =", S.N(minus, 20))

print("R29I PASS")
print("TWO OFF-CENTER OBSTRUCTIONS: EXACTLY CERTIFIED")
print("GENERAL FAMILY CLOSURE: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
