import itertools
import sympy as sp

print("Checkpoint 016: Full-matrix robustness")

a,b,c,d = sp.symbols("a b c d", real=True)

A = sp.Matrix([
    [-1+a, 10+b],
    [c, -2+d]
])

trace = sp.trace(A)
det = A.det()

assert sp.expand(trace-(-3+a+d)) == 0
assert sp.expand(
    det-(2-2*a-d+a*d-10*c-b*c)
) == 0

v = sp.Matrix([1,1])/sp.sqrt(2)
growth = sp.expand(2*(v.T*A*v)[0])

assert sp.simplify(
    growth-(7+a+b+c+d)
) == 0

print("PASS: exact symbolic identities")

eps = sp.Rational(1,20)

# Conservative bounds valid for every
# perturbation in the entire box.
trace_upper = -3+2*eps

det_lower = (
    2-2*eps-eps-eps**2
    -10*eps-eps**2
)

growth_lower = 7-4*eps

assert trace_upper < 0
assert det_lower > 0
assert growth_lower > 0

print("PASS: global analytical bounds")
print("trace_upper =", trace_upper)
print("det_lower =", det_lower)
print("growth_lower =", growth_lower)

# Independent finite grid verification.
values = [-eps,sp.Integer(0),eps]
tested = 0

for entries in itertools.product(values,repeat=4):
    M = A.subs(dict(zip((a,b,c,d),entries)))

    assert sp.trace(M) < 0
    assert M.det() > 0

    rate = sp.simplify(2*(v.T*M*v)[0])
    assert rate > 0

    tested += 1

assert tested == 81

print("grid_matrices_tested =", tested)
print("PASS: full-matrix robustness verified")
