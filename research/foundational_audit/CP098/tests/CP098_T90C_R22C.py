import sympy as s

print("=== CP098-T90C-R22C ===")

eta = s.diag(-1, 1, 1, 1)

def dot(x, y):
    return (x.T * eta * y)[0]

u = s.Matrix([5, 3, 0, 0])
q = s.Matrix([1, 2, 2, 1])
j = s.Matrix([4, 1, 3, 6])

A = s.Integer(3)
C = s.Integer(2)

p0 = s.Matrix([
    s.Rational(1797334885453, 504122720256),
    -s.Rational(7603294311121, 1512368160768),
    s.Rational(7357855225, 7876917504),
    s.Rational(20080631585, 15753835008)
])

v = s.sqrt(-dot(u, u))
T = u/v
w = q + dot(T, q)*T
w2 = s.simplify(dot(w, w))
k = s.sqrt(w2)/v**2
N = w/s.sqrt(w2)

projector = (
    s.eye(4)
    + T*(T.T*eta)
    - N*(N.T*eta)
)

R = projector*j

P = -2*C*eta*R/(v*w2)

Lgeom = (
    -A*w2/v**3
    -C*dot(R, R)/(v*w2)
)

# R20: p1.T=0 and p1.N=-2*A*k/v.
# q lies in the span of T and N.
p1q = (-2*A*k/v)*dot(N, q)

p0u = s.simplify((p0.T*u)[0])
Pj = s.simplify((P.T*j)[0])

H = s.simplify(p0u+p1q+Pj-Lgeom)

print("p0.u =", p0u)
print("p1.q =", s.simplify(p1q))
print("P.j =", Pj)
print("Lgeom =", s.simplify(Lgeom))
print("Hgeom =", H)

assert H == 0

print()
print("R22C PASS")
print("GEOMETRIC HAMILTONIAN IDENTITY: VERIFIED AT ONE JET")
print("INDEPENDENT FULL p0 DERIVATION: STILL OPEN")
print("FULL ACTION AND DIRAC SYSTEM: OPEN")
