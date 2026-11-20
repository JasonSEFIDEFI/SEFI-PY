import sympy as S
from time import perf_counter

t0 = perf_counter()
print("=== CP098-T90C-R27S ===", flush=True)

eta = S.diag(-1, 1, 1, 1)
u = S.Matrix([2, 0, 0, 0])
q = S.Matrix([0, 1, 0, 0])
P = S.Matrix([0, 0, 1, 0])
p1 = S.zeros(4, 1)

def dot(x, y):
    return (x.T * eta * y)[0]

s = -dot(u, u)
v = S.sqrt(s)
d = dot(u, q)
W = dot(q, q) + d**2/s
P2 = dot(P, P)

IA = -S.Rational(3, 2)
IW = S.Integer(1)
hA = hS = hW = S.Integer(1)
kA = kS = kW = S.Integer(1)
m = S.Integer(1)
pA = pS = pW = S.Integer(0)

IS = (
    3*kA*v**4/(4*W)
    +28*P2*v**4/(961*IW)
)

A = (3+2*IA)*IS
C = (2+IA+S.Rational(3,2)*IA**2)*IW

V = (kA*IA**2+kS*IS**2+kW*IW**2)/2
T = (pA**2/hA+pS**2/hS+pW**2/hW)/2
M = m+V+T
J = 2*IS*pA/hA+(3+2*IA)*pS/hS
b = 2*A*W/v**3

# Set chi=psi=theta=B=0 by construction.
p0 = S.Matrix([
    S.factor(v*W*P2/(4*C)/2-v*M/2),
    0, 0, 0
])

phi1 = (u.T*P)[0]
phi2 = (q.T*P)[0]
chi = (u.T*p1)[0]
psi = (q.T*p1)[0]+b

theta = (
    -(p0.T*u)[0]
    +(p1.T*q)[0]
    +v*W*P2/(4*C)
    +3*A*W/v**3-v*M
)

B = (
    -(p0.T*q)[0]
    -v*W*dot(p1,P)/(2*C)
    +d*W*P2/(4*C*v)
    +5*A*d*W/v**5
    +d*M/v
    -2*A*W*(phi2+d*phi1/s)/(C*v**2)
    +2*W*J/v**2
)

Q = (
    S.Rational(3,2)*v*kA
    +56*v*W*P2/(961*IW)
    -2*IS*W/v**3
)

print("STAGE 1: regular domain")
for name, value in [
    ("s", s), ("v", v), ("W", W),
    ("C", C), ("P²", P2), ("IS", IS)
]:
    print(name, "=", value)
assert s > 0 and W > 0 and C != 0
print("STAGE 1 PASS")

print("STAGE 2: exact constraints")
for name, value in [
    ("phi1", phi1), ("phi2", phi2),
    ("chi", chi), ("psi", psi),
    ("theta", theta), ("B", B),
    ("A", A), ("J", J), ("Q", Q),
    ("pA", pA)
]:
    residual = S.cancel(value)
    if residual != 0:
        raise AssertionError(f"{name} = {residual}")
    print("PASS:", name)

assert IS != 0
print("STAGE 2 PASS")

print()
print("R27S COMPLETE")
print("EXACT 4D BRANCH WITNESS: CONSTRUCTED")
print("D AT WITNESS: NOT YET EVALUATED")
print("BRANCH PRESERVATION: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(perf_counter()-t0, 3))
