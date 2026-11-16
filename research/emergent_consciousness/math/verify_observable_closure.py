import sympy as sp

a, b, u, v1, v2 = sp.symbols("a b u v1 v2")

next_1 = a*u + b*v1
next_2 = a*u + b*v2

difference = sp.factor(next_1 - next_2)

assert sp.simplify(difference - b*(v1-v2)) == 0
assert sp.simplify(difference.subs(b, 0)) == 0

print("Observable closure verification")
print("Difference:", difference)
print("b = 0: observable evolution is closed")
print("b != 0: distinct hidden states can yield distinct futures")
print("PASS: symbolic conditions verified")
