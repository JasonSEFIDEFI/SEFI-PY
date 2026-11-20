import sympy as s
from time import monotonic

print("=== CP098-T90C-R24B ===", flush=True)
start = monotonic()

# Univariate second-order Taylor algebra.
# A polynomial represents f0 + f1*t + f2*t**2.
# Thus the second derivative at zero is 2*f2.

class T2:
    def __init__(self, a=0, b=0, c=0):
        self.x = (s.sympify(a), s.sympify(b), s.sympify(c))

    def __add__(self, other):
        other = cv(other)
        return T2(*(self.x[i]+other.x[i] for i in range(3)))

    __radd__ = __add__

    def __neg__(self):
        return T2(*(-x for x in self.x))

    def __sub__(self, other):
        return self + (-cv(other))

    def __rsub__(self, other):
        return cv(other) - self

    def __mul__(self, other):
        other = cv(other)
        a,b,c = self.x
        d,e,f = other.x
        return T2(a*d, a*e+b*d, a*f+b*e+c*d)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        p = s.sympify(exponent)
        a,b,c = self.x

        if a == 0:
            raise ZeroDivisionError("STOP: zero Taylor base")

        return T2(
            a**p,
            p*a**(p-1)*b,
            p*a**(p-1)*c
            + p*(p-1)*a**(p-2)*b*b/2
        )

    def __truediv__(self, other):
        return self * (cv(other)**(-1))

    def __rtruediv__(self, other):
        return cv(other) * (self**(-1))

def cv(x):
    return x if isinstance(x,T2) else T2(x)

def add(x,y):
    return [a+b for a,b in zip(x,y)]

def scale(c,x):
    return [c*a for a in x]

def dot(x,y):
    return -x[0]*y[0] + sum(
        (x[i]*y[i] for i in range(1,4)), T2()
    )

def evaluate(u0,q0,j0,a0,b0,C0,C1,C2):
    u = [
        T2(u0[i],q0[i],s.Rational(j0[i],2))
        for i in range(4)
    ]
    q = [
        T2(q0[i],j0[i],s.Rational(a0[i],2))
        for i in range(4)
    ]
    j = [
        T2(j0[i],a0[i],s.Rational(b0[i],2))
        for i in range(4)
    ]

    C = T2(C0,C1,s.Rational(C2,2))

    v2 = -dot(u,u)
    v = v2**s.Rational(1,2)

    # T = u/v; w is orthogonal to T.
    T = scale(1/v,u)
    w = add(q,scale(dot(T,q),T))
    w2 = dot(w,w)

    # R = transverse projection of j.
    R = add(
        add(j,scale(dot(T,j),T)),
        scale(-dot(w,j)/w2,w)
    )

    F = scale(1/(v*w2),R)

    # eta lowers the index: (-,+,+,+).
    F[0] = -F[0]

    # Product-rule decomposition, computed exactly.
    terms = []
    for f in F:
        f0,f1,f2 = f.x

        c0 = s.sympify(C0)
        c1 = s.sympify(C1)
        c2 = s.sympify(C2)

        result = -2*(
            c2*f0 + 2*c1*f1 + 2*c0*f2
        )
        terms.append(s.factor(result))

    return s.Matrix(terms)

cases = [
    (
        "JET 1",
        [5,3,0,0],
        [1,2,2,1],
        [4,1,3,6],
        [2,-1,4,3],
        [1,3,-2,5]
    ),
    (
        "JET 2",
        [13,3,4,0],
        [2,3,1,4],
        [1,4,2,5],
        [3,-2,1,4],
        [2,1,-3,5]
    )
]

for index,(name,u,q,j,a,b) in enumerate(cases):
    print()
    print(name, "START", flush=True)

    result = evaluate(u,q,j,a,b,2,1,3)

    print("P SECOND DERIVATIVE =", result, flush=True)

    if index == 0:
        # Exact R24A/R22B value, obtained by a
        # different symbolic differentiation method.
        expected0 = s.Rational(3250183453,492307344)

        residual = s.simplify(result[0]-expected0)

        print("R24A COMPONENT-0 RESIDUAL =", residual)
        assert residual == 0

        expected_vector = s.Matrix([
            s.Rational(3250183453,492307344),
            -s.Rational(18580716433,1476922032),
            -s.Rational(248659083,54700816),
            s.Rational(110869279,492307344)
        ])

        residual_vector = (
            result-expected_vector
        ).applyfunc(s.simplify)

        print("R22B VECTOR RESIDUAL =", residual_vector)
        assert residual_vector == s.zeros(4,1)

    print(name, "PASS", flush=True)
    print("ELAPSED SECONDS =", round(monotonic()-start,2),
          flush=True)

print()
print("R24B COMPLETE")
print("EXACT P SECOND DERIVATIVE: TWO JETS")
print("FIRST-JET COMPONENT COMPARISON: PASS")
print("FULL p0 COMPONENT COMPARISON: STILL OPEN")
print("FULL ACTION AND DIRAC SYSTEM: OPEN")
