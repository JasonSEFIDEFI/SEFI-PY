import sympy as s

print("=== CP098-T90C-R22B ===", flush=True)

# Exact arithmetic in Q[t,e]/(t^3,e^2).
# Coefficients are indexed by (power of t, power of e).
# No floating-point approximation is used.

class Jet:
    def __init__(self, terms=None):
        self.c = {
            key: s.sympify(value)
            for key, value in (terms or {}).items()
            if value != 0
        }

    @staticmethod
    def const(x):
        return Jet({(0, 0): x})

    def __add__(self, other):
        other = asjet(other)
        out = dict(self.c)
        for key, value in other.c.items():
            out[key] = out.get(key, 0) + value
        return Jet(out)

    __radd__ = __add__

    def __neg__(self):
        return Jet({k: -v for k, v in self.c.items()})

    def __sub__(self, other):
        return self + (-asjet(other))

    def __rsub__(self, other):
        return asjet(other) - self

    def __mul__(self, other):
        other = asjet(other)
        out = {}
        for (i, j), x in self.c.items():
            for (k, l), y in other.c.items():
                if i + k <= 2 and j + l <= 1:
                    key = (i + k, j + l)
                    out[key] = out.get(key, 0) + x*y
        return Jet(out)

    __rmul__ = __mul__

    def __pow__(self, power):
        power = s.sympify(power)
        base = self.c.get((0, 0), s.S.Zero)
        if base == 0:
            raise ZeroDivisionError(
                "STOP: nonregular jet or zero expansion base"
            )

        nil = self - base
        z = nil * (1/base)

        result = Jet.const(1)
        term = Jet.const(1)

        # Nilpotent fourth power vanishes in this algebra.
        for n in range(1, 4):
            term = term*z
            result = result + s.binomial(power, n)*term

        return (base**power)*result

    def __truediv__(self, other):
        return self * (asjet(other)**(-1))

    def __rtruediv__(self, other):
        return asjet(other) * (self**(-1))

    def coeff(self, tpower, epower):
        return s.simplify(
            self.c.get((tpower, epower), s.S.Zero)
        )

def asjet(x):
    return x if isinstance(x, Jet) else Jet.const(x)

def vec(values):
    return [asjet(x) for x in values]

def vadd(x, y):
    return [a+b for a, b in zip(x, y)]

def scale(c, x):
    return [c*a for a in x]

def mdot(x, y):
    return -x[0]*y[0] + sum(
        (x[i]*y[i] for i in range(1, 4)),
        Jet.const(0)
    )

def lagrangian(u, q, j, A, C):
    v2 = -mdot(u, u)
    v = v2**s.Rational(1, 2)

    T = scale(1/v, u)
    w = vadd(q, scale(mdot(T, q), T))
    w2 = mdot(w, w)

    R = vadd(
        vadd(j, scale(mdot(T, j), T)),
        scale(-mdot(w, j)/w2, w)
    )

    R2 = mdot(R, R)

    return -A*w2/(v**3) - C*R2/(v*w2)

u0 = [5, 3, 0, 0]
q0 = [1, 2, 2, 1]
j0 = [4, 1, 3, 6]
a0 = [2, -1, 4, 3]
b0 = [1, 3, -2, 5]

A0, A1 = s.Integer(3), s.Integer(2)
C0, C1, C2 = s.Integer(2), s.Integer(1), s.Integer(3)

def trajectory(b):
    u = [
        Jet({(0,0):u0[i], (1,0):q0[i],
             (2,0):s.Rational(j0[i],2)})
        for i in range(4)
    ]
    q = [
        Jet({(0,0):q0[i], (1,0):j0[i],
             (2,0):s.Rational(a0[i],2)})
        for i in range(4)
    ]
    j = [
        Jet({(0,0):j0[i], (1,0):a0[i],
             (2,0):s.Rational(b[i],2)})
        for i in range(4)
    ]
    A = Jet({(0,0):A0, (1,0):A1})
    C = Jet({(0,0):C0, (1,0):C1,
             (2,0):C2/2})
    return u, q, j, A, C

def perturb(x, i):
    result = list(x)
    result[i] = result[i] + Jet({(0,1):1})
    return result

def momentum_parts(b, c2=C2, progress=False):
    u, q, j, A, C = trajectory(b)
    C = Jet({(0,0):C0, (1,0):C1,
             (2,0):c2/2})

    if progress:
        print("STAGE 1: partial-u gradient", flush=True)

    gu = s.Matrix([
        lagrangian(perturb(u,i),q,j,A,C).coeff(0,1)
        for i in range(4)
    ])

    if progress:
        print("STAGE 1 COMPLETE", flush=True)
        print("STAGE 2: total q-gradient derivative", flush=True)

    gqdot = s.Matrix([
        lagrangian(u,perturb(q,i),j,A,C).coeff(1,1)
        for i in range(4)
    ])

    if progress:
        print("STAGE 2 COMPLETE", flush=True)
        print("STAGE 3: second total P derivative", flush=True)

    Pddot = s.Matrix([
        2*lagrangian(u,q,perturb(j,i),A,C).coeff(2,1)
        for i in range(4)
    ])

    if progress:
        print("STAGE 3 COMPLETE", flush=True)

    p0 = (gu-gqdot+Pddot).applyfunc(s.simplify)
    return gu, gqdot, Pddot, p0

gu, gqdot, Pddot, p0 = momentum_parts(
    b0, progress=True
)

print()
print("PARTIAL-u GRADIENT =", gu)
print("q-GRADIENT DERIVATIVE =", gqdot)
print("SECOND P DERIVATIVE =", Pddot)
print("GEOMETRIC p0 =", p0, flush=True)

# Independent structural check:
# Changing only the fifth derivative must change p0 by W*delta_b.

print()
print("STAGE 4: independent fifth-derivative response", flush=True)

zero_b = [0,0,0,0]
p0_zero_b = momentum_parts(zero_b)[3]

eta = s.diag(-1,1,1,1)
U = s.Matrix(u0)
Q = s.Matrix(q0)
B = s.Matrix(b0)

def dot(x,y):
    return (x.T*eta*y)[0]

v = s.sqrt(-dot(U,U))
T = U/v
w = Q + dot(T,Q)*T
w2 = s.simplify(dot(w,w))
N = w/s.sqrt(w2)

projector = (
    s.eye(4)
    + T*(T.T*eta)
    - N*(N.T*eta)
)

W = -2*C0*eta*projector/(v*w2)

actual = (p0-p0_zero_b).applyfunc(s.simplify)
expected = (W*B).applyfunc(s.simplify)

assert (actual-expected).applyfunc(
    s.simplify
) == s.zeros(4,1)

print("FIFTH-DERIVATIVE RESIDUAL =",
      (actual-expected).applyfunc(s.simplify))
print("PASS", flush=True)

# Independent coupling-derivative check:
# C'' enters p0 through the second derivative of P.
print("STAGE 5: C-second-derivative response", flush=True)

p0_zero_c2 = momentum_parts(b0, c2=0)[3]

J = s.Matrix(j0)
R = projector*J
P0 = -2*C0*eta*R/(v*w2)

actual_c2 = (p0-p0_zero_c2).applyfunc(s.simplify)
expected_c2 = (C2/C0*P0).applyfunc(s.simplify)

assert (actual_c2-expected_c2).applyfunc(
    s.simplify
) == s.zeros(4,1)

print("C-SECOND-DERIVATIVE RESIDUAL =",
      (actual_c2-expected_c2).applyfunc(s.simplify))
print("PASS", flush=True)

print()
print("R22B COMPLETE")
print("GEOMETRIC p0: EXACTLY EVALUATED AT ONE JET")
print("FIFTH-DERIVATIVE RESPONSE: VERIFIED")
print("C-SECOND-DERIVATIVE RESPONSE: VERIFIED")
print("SECOND INDEPENDENT JET: PENDING")
print("FULL ACTION MOMENTUM: OPEN")
print("HAMILTONIAN RECONSTRUCTION: OPEN")
print("FULL DIRAC SYSTEM: OPEN")
print("DEFI TORSION IDENTIFICATION: UNRESOLVED")
