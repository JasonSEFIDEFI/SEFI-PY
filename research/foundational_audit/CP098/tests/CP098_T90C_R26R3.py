import sympy as S
from time import monotonic

t0 = monotonic()

print("=== CP098-T90C-R26R3 ===", flush=True)

def check(label, residual):
    result = S.factor(residual)
    if result != 0:
        print("FAIL:", label, flush=True)
        print("RESIDUAL:", result, flush=True)
        raise AssertionError(label)
    print("PASS:", label, flush=True)

# Formal geometric invariants.
# s>0, W>0, C!=0, h_a>0.
#
# v=sqrt(s), d=u.eta.q.
#
# A and C depend only on identity variables.

v, W, d, s = S.symbols(
    "v W d s", positive=True
)

A, C, M = S.symbols(
    "A C M", nonzero=True
)

J, K, P2, p1P = S.symbols(
    "J K P2 p1P"
)

phi1, phi2 = S.symbols(
    "phi1 phi2"
)

# Identity derivatives and momenta.
h = S.symbols("hA hS hW", positive=True)
pI = S.symbols("pIA pIS pIW")
dA = S.symbols("dAA dAS dAW")

# Canonical Hamiltonian pieces:
#
# Hc = K + J + D + F + G
#
# K = p0.u
# J = p1.q
# D = -vW P²/(4C)
# F = AW/v³
# G = vM

D = -v*W*P2/(4*C)
F = A*W/v**3
G = v*M

E = 2*F

print("STAGE 1: geometric differential rules", flush=True)

# Directional derivatives along q in u-space.
#
# q.du(v) = -d/v
# q.du(W) = 2dW/s
#
# s=v².

check(
    "s=v² substitution",
    (2*d*W/s).subs(s, v**2)
    - 2*d*W/v**2
)

# q.du(vW) = dW/v

qdu_vW = (
    (-d/v)*W
    + v*(2*d*W/v**2)
)

check(
    "q.du(vW)=dW/v",
    qdu_vW-d*W/v
)

# q.du(AW/v³)=5AdW/v5

qdu_F = A*(
    (2*d*W/v**2)/v**3
    - 3*W*(-d/v)/v**4
)

check(
    "q.du(F)=5AdW/v5",
    qdu_F-5*A*d*W/v**5
)

# q.du(E)=2*q.du(F)

check(
    "q.du(E)=10AdW/v5",
    2*qdu_F-10*A*d*W/v**5
)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: negative u/p1 contribution", flush=True)

# -q.du(Hc)
#
# q.du(K) = p0.q, denoted Kq.
# q.du(J) = 0.
#
# q.du(D) = -dW P²/(4Cv).
# q.du(F) = 5AdW/v5.
# q.du(G) = -dM/v.

Kq = S.symbols("p0q")

qdu_D = -d*W*P2/(4*C*v)
qdu_G = -d*M/v

negative_direct = (
    -Kq
    -qdu_D
    -qdu_F
    -qdu_G
)

negative_expected = (
    -Kq
    +d*W*P2/(4*C*v)
    -5*A*d*W/v**5
    +d*M/v
)

check(
    "negative u/p1 term",
    negative_direct-negative_expected
)

print("STAGE 2 PASS", flush=True)

print("STAGE 3: q/P contribution", flush=True)

# dHc/dP = -vW/(2C) eta.P
#
# dpsi/dq = p1 + dE/dq.
#
# dE/dq = 2A/v³ * dW/dq.
#
# dW/dq = 2 eta.q + 2d/s eta.u.
#
# Thus:
# dE/dq . eta.P
# = 4A/v³*(phi2+d phi1/s).

qP_direct = (
    -v*W*p1P/(2*C)
    -2*A*W/(C*v**2)
       *(phi2+d*phi1/v**2)
)

qP_expected = (
    -v*W*p1P/(2*C)
    -2*A*W/(C*v**2)
       *(phi2+d*phi1/s)
)

check(
    "q/P contribution",
    (qP_direct-qP_expected).subs(s,v**2)
)

print("STAGE 3 PASS", flush=True)

print("STAGE 4: identity-sector contribution", flush=True)

# dpsi/dI_a = 2W/v³ * dA/dI_a
#
# dHc/dpI_a = v*pI_a/h_a

identity_direct = sum(
    (2*W/v**3)*dA[i]
    *(v*pI[i]/h[i])
    for i in range(3)
)

identity_expected = (
    2*W/v**2
    *sum(dA[i]*pI[i]/h[i] for i in range(3))
)

check(
    "identity-sector contribution",
    identity_direct-identity_expected
)

print("STAGE 4 PASS", flush=True)

print("STAGE 5: assemble general B", flush=True)

positive_term = 10*A*d*W/v**5

B_structured = (
    positive_term
    +negative_direct
    +qP_direct
    +identity_direct
)

B_expected = (
    -Kq
    -v*W*p1P/(2*C)
    +d*W*P2/(4*C*v)
    +5*A*d*W/v**5
    +d*M/v
    -2*A*W/(C*v**2)
       *(phi2+d*phi1/v**2)
    +identity_expected
)

check(
    "GENERAL B FORMULA",
    B_structured-B_expected
)

print("STAGE 5 PASS", flush=True)

print()
print("R26R3 COMPLETE")
print("STRUCTURED GENERAL B DERIVATION: PASS")
print("IDENTITY COUPLING: INCLUDED")
print("NO NUMERICAL WITNESSES USED")
print("FULL DIRECT SYMBOLIC DIFFERENTIATION: OPEN")
print("CONSTRAINT REGULARITY: OPEN")
print("A-ZERO BRANCH: OPEN")
print("GAUGE GENERATORS: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(monotonic()-t0, 3)
)
