import sympy as sp

a, b, c, d = sp.symbols("a b c d")
u0, v0 = sp.symbols("u0 v0")

A = sp.Matrix([[a, b], [c, d]])
x0 = sp.Matrix([u0, v0])

x1 = A * x0
x2 = A * x1

u1 = x1[0]
u2 = x2[0]

predicted = (a + d)*u1 + (b*c - a*d)*u0
residual = sp.simplify(u2 - predicted)

print("Exact recurrence verification")
print("Residual:", residual)

assert residual == 0
print("PASS: symbolic identity verified")
