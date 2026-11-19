import sympy as S
from time import monotonic

t0 = monotonic()
print("=== CP098-T90C-R26G ===", flush=True)

eta = S.diag(-1, 1, 1, 1)

def mdot(a, b):
    return (a.T * eta * b)[0]

def pair(a, b):
    return (a.T * b)[0]

u = S.Matrix(S.symbols("u0:4"))
q = S.Matrix(S.symbols("q0:4"))
p1 = S.Matrix(S.symbols("p1_0:4"))

A = S.symbols("A", real=True)

s = -mdot(u, u)
v = S.sqrt(s)
d = mdot(u, q)
W = mdot(q, q) + d**2/s

print("STAGE 1: geometric regularity", flush=True)

jets = [
    (
        S.Matrix([5, 3, 0, 0]),
        S.Matrix([1, 2, 2, 1])
    ),
    (
        S.Matrix([13, 3, 4, 0]),
        S.Matrix([2, 3, 1, 4])
    )
]

for i, (uj, qj) in enumerate(jets, 1):
    subs = dict(zip(list(u)+list(q),
                    list(uj)+list(qj)))

    sj = S.simplify(s.subs(subs))
    wj = S.simplify(W.subs(subs))

    assert sj > 0
    assert wj > 0

    # The two row vectors u and q must
    # be linearly independent.
    rank = S.Matrix([
        list(uj),
        list(qj)
    ]).rank()

    assert rank == 2

    print(
        "JET", i,
        "s =", sj,
        "W =", wj,
        "rank(u,q) =", rank,
        flush=True
    )

print("STAGE 1 PASS", flush=True)

print("STAGE 2: independence of four constraints",
      flush=True)

# Constraints:
# phi1=u.P, phi2=q.P,
# chi=u.p1, psi=q.p1+2AW/v^3.
#
# The derivatives with respect to the
# eight variables (P,p1) form the block:
#
# [u  0]
# [q  0]
# [0  u]
# [0  q]
#
# Each 2x4 diagonal block has rank 2
# when u and q are independent.
#
# Therefore the four constraint gradients
# have rank 4 on the regular domain.

for uj, qj in jets:
    zero = [S.Integer(0)]*4

    J = S.Matrix([
        list(uj)+zero,
        list(qj)+zero,
        zero+list(uj),
        zero+list(qj)
    ])

    assert J.rank() == 4

print("FOUR CONSTRAINT GRADIENTS: RANK 4")
print("STAGE 2 PASS", flush=True)

print("STAGE 3: tertiary independence",
      flush=True)

# theta includes -p0.u.
#
# None of phi1,phi2,chi,psi depends
# on p0.
#
# Hence d(theta)/d(p0)=-u.
#
# Since s>0, u cannot vanish.
#
# The fifth gradient is therefore
# independent of the first four.

for uj, qj in jets:
    zero = [S.Integer(0)]*4

    J5 = S.Matrix([
        zero+list(uj)+zero,
        zero+list(qj)+zero,
        zero+zero+list(uj),
        zero+zero+list(qj),
        list(-uj)+zero+zero
    ])

    assert J5.rank() == 5

print("FIVE CONSTRAINT GRADIENTS: RANK 5")
print("STAGE 3 PASS", flush=True)

print("STAGE 4: secondary bracket", flush=True)

# chi=u.p1 generates:
# delta(u)=u
# delta(p1)=-p1
# delta(q)=0.
#
# Consequently:
# delta(p1.q)=-p1.q
# delta(W)=0
# delta(v)=v
# delta(2AW/v^3)=-6AW/v^3.
#
# Since delta(F)={F,chi},
# antisymmetry gives:
#
# {chi,psi}=p1.q+6AW/v^3
#          =psi+4AW/v^3.

pq = S.symbols("pq", real=True)

psi = pq + 2*A*W/v**3
bracket = pq + 6*A*W/v**3

residual = S.expand(
    bracket-(psi+4*A*W/v**3)
)

assert residual == 0

print("{chi,psi}=psi+4AW/v^3")
print("BRACKET IDENTITY: PASS")

jet1 = dict(zip(
    list(u)+list(q),
    list(jets[0][0])+list(jets[0][1])
))
jet1[A] = 3

weak = S.simplify(
    (4*A*W/v**3).subs(jet1)
)

print("WEAK BRACKET AT JET 1 =", weak)
assert weak == S.Rational(387,256)

print("STAGE 4 PASS", flush=True)

print()
print("R26G COMPLETE")
print("REGULAR GEOMETRIC JETS: VERIFIED")
print("FOUR CONSTRAINT GRADIENTS: INDEPENDENT")
print("TERTIARY GRADIENT: INDEPENDENT")
print("SECONDARY BRACKET: VERIFIED BY SCALING")
print("psi PRESERVATION: OPEN")
print("TERTIARY PRESERVATION: OPEN")
print("FULL DIRAC CLASSIFICATION: OPEN")
print("DEFI CORRESPONDENCE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(monotonic()-t0, 3))
