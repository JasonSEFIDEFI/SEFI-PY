import sympy as S
from time import monotonic

t0 = monotonic()
print("=== CP098-T90C-R26E ===", flush=True)

eta = S.diag(-1, 1, 1, 1)

u = S.Matrix(S.symbols("u0:4"))
q = S.Matrix(S.symbols("q0:4"))
P = S.Matrix(S.symbols("P0:4"))
p1 = S.Matrix(S.symbols("p1_0:4"))

A, C = S.symbols("A C", nonzero=True)

def mdot(a, b):
    return (a.T * eta * b)[0]

def pair(a, b):
    return (a.T * b)[0]

s = -mdot(u, u)
v = S.sqrt(s)
d = mdot(u, q)
W = mdot(q, q) + d**2/s
P2 = mdot(P, P)

phi1 = pair(u, P)
phi2 = pair(q, P)
chi = pair(u, p1)
psi = pair(q, p1) + 2*A*W/v**3

def gradient(f, variables):
    return S.Matrix([
        S.diff(f, x) for x in variables
    ])

def PB(f, g):
    result = S.Integer(0)

    for x, p in (
        list(zip(u, p1)) + list(zip(q, P))
    ):
        result += (
            S.diff(f, x)*S.diff(g, p)
            - S.diff(f, p)*S.diff(g, x)
        )

    return result

def check(name, actual, expected):
    residual = S.simplify(actual-expected)
    print(name, "RESIDUAL =", residual, flush=True)
    assert residual == 0, name

print("STAGE 1: W derivative identities", flush=True)

gradW = gradient(W, q)

check("u.grad_q(W)", pair(u, gradW), 0)
check("q.grad_q(W)", pair(q, gradW), 2*W)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: independent primary preservation",
      flush=True)

# Terms of H independent of q and P do not
# contribute to either primary bracket.
#
# H = p0.u + p1.q
#     - v*W*P2/(4*C)
#     + A*W/v**3
#     + terms depending on u and identity data.

Hrelevant = (
    pair(p1, q)
    - v*W*P2/(4*C)
    + A*W/v**3
)

check(
    "{phi1,H}",
    PB(phi1, Hrelevant),
    phi2-chi
)

check(
    "{phi2,H}",
    PB(phi2, Hrelevant),
    -psi
)

print("STAGE 2 PASS", flush=True)

print("STAGE 3: primary-secondary brackets",
      flush=True)

check(
    "{phi1,chi}",
    PB(phi1, chi),
    phi1
)

check(
    "{phi2,chi}",
    PB(phi2, chi),
    0
)

check(
    "{phi1,psi}",
    PB(phi1, psi),
    phi2-chi
)

check(
    "{phi2,psi}",
    PB(phi2, psi),
    -psi-2*A*W/v**3
)

print("STAGE 3 PASS", flush=True)

print("STAGE 4: nonzero weak bracket witness",
      flush=True)

# On the full candidate constraint surface:
# phi1=phi2=chi=psi=0.
#
# The final primary-secondary bracket reduces to:
# {phi2,psi} = -2*A*W/v**3.

weak_value = -2*A*W/v**3

jet = {
    u[0]: 5, u[1]: 3,
    u[2]: 0, u[3]: 0,
    q[0]: 1, q[1]: 2,
    q[2]: 2, q[3]: 1,
    A: 3, C: 2
}

wv = S.simplify(W.subs(jet))
vv = S.simplify(v.subs(jet))
bv = S.simplify(weak_value.subs(jet))

print("JET v =", vv)
print("JET W =", wv)
print("WEAK BRACKET =", bv)

assert vv == 4
assert wv == S.Rational(129, 16)
assert bv == -S.Rational(387, 1024)
assert bv != 0

print("STAGE 4 PASS", flush=True)

print()
print("R26E COMPLETE")
print("W DIFFERENTIAL IDENTITIES: VERIFIED")
print("PRIMARY PRESERVATION: INDEPENDENTLY VERIFIED")
print("PRIMARY-SECONDARY BRACKETS: VERIFIED")
print("NONZERO WEAK BRACKET: VERIFIED")
print("ALL-FIRST-CLASS CLASSIFICATION: EXCLUDED")
print("SECONDARY PRESERVATION: OPEN")
print("FULL DIRAC CLASSIFICATION: OPEN")
print("PHYSICAL DEGREES OF FREEDOM: NOT COUNTED")
print("DEFI CORRESPONDENCE: OPEN")
print("ELAPSED SECONDS =",
      round(monotonic()-t0, 3))
