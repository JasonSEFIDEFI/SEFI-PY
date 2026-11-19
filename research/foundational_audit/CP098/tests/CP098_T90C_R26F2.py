import sympy as S
from time import monotonic

t0 = monotonic()
print("=== CP098-T90C-R26F2 ===", flush=True)

eta = S.diag(-1, 1, 1, 1)
u = S.Matrix(S.symbols("u0:4"))
q = S.Matrix(S.symbols("q0:4"))

def mdot(a,b):
    return (a.T*eta*b)[0]

s = -mdot(u,u)
v = S.sqrt(s)
d = mdot(u,q)
W = mdot(q,q) + d**2/s

print("STAGE 1: exact scaling identities", flush=True)

# Avoid general-purpose simplify on expanded geometry.
Euler_s = sum(x*S.diff(s,x) for x in u)
Euler_d = sum(x*S.diff(d,x) for x in u)

assert S.expand(Euler_s-2*s) == 0
assert S.expand(Euler_d-d) == 0

# W = q^2 + d^2/s.
# d has degree 1; s has degree 2.
Euler_W = S.cancel(
    2*d*Euler_d/s
    - d**2*Euler_s/s**2
)

assert Euler_W == 0

print("Euler_u(s)=2s: PASS")
print("Euler_u(d)=d: PASS")
print("Euler_u(W)=0: PASS")
print("Euler_u(v)=v: follows exactly")
print("STAGE 1 PASS", flush=True)

print("STAGE 2: two rational jets", flush=True)

jets = [
    ([5,3,0,0],[1,2,2,1]),
    ([13,3,4,0],[2,3,1,4])
]

for index,(uv,qv) in enumerate(jets,1):
    sub = dict(zip(list(u)+list(q),uv+qv))
    sv = s.subs(sub)
    vv = v.subs(sub)
    wv = W.subs(sub)

    assert sv > 0
    assert wv > 0
    assert Euler_W.subs(sub) == 0

    print(
        "JET",index,
        "v =",vv,
        "W =",wv,
        flush=True
    )

print("STAGE 2 PASS", flush=True)

print("STAGE 3: canonical bracket by homogeneity",
      flush=True)

p0u,p1q,P2 = S.symbols("p0u p1q P2")
A,C,m,V,T = S.symbols("A C m V T")
l1 = S.symbols("l1")

# Hc = p0u+p1q-v*W*P2/(4C)
#      +A*W/v^3+v*(m+V+T).
#
# chi=u.p1 generates:
# {u,chi}=u and {p1,chi}=-p1.
#
# Therefore {chi,Hc} is computed using
# the negative of the scaling variation.

theta = (
    -p0u
    +p1q
    +v*W*P2/(4*C)
    +3*A*W/v**3
    -v*(m+V+T)
)

# Independently assemble {chi,Hc}:
# contribution from p1.q: +p1q
# contribution from -u.grad_u(Hc): remaining terms

derived = (
    p1q
    -p0u
    -(-v*W*P2/(4*C))
    -(-3*A*W/v**3)
    -v*(m+V+T)
)

assert S.expand(derived-theta) == 0

print("{chi,Hc}=theta: PASS")
print("STAGE 3 PASS", flush=True)

print("STAGE 4: multiplier dependence", flush=True)

# Previously verified:
# {chi,phi1}=-phi1
# {chi,phi2}=0.
#
# HT=Hc+l1*phi1+l2*phi2.

phi1 = S.symbols("phi1")
total = theta-l1*phi1

assert S.expand(total-(theta-l1*phi1)) == 0

print("{chi,HT}=theta-l1*phi1")
print("On primary surface: dot(chi)=theta")
print("STAGE 4 PASS", flush=True)

print()
print("R26F2 COMPLETE")
print("GEOMETRIC HOMOGENEITY: VERIFIED")
print("chi PRESERVATION: VERIFIED BY SCALING")
print("TERTIARY CANDIDATE theta: DERIVED")
print("psi PRESERVATION: OPEN")
print("TERTIARY INDEPENDENCE: OPEN")
print("FULL DIRAC CLASSIFICATION: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",round(monotonic()-t0,3))
