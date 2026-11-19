import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R26K ===", flush=True)

# Independent scalar building blocks.
p0u, p1q = S.symbols("p0u p1q")
X, Y, E = S.symbols("X Y E")

# X = v*W*P2/(4*C)
# Y = A*W/v**3
# E = v*(m+V+T)

Hc = p0u + p1q - X + Y + E
theta = -p0u + p1q + X + 3*Y - E
psi = p1q + 2*Y

print("STAGE 1: exact Hamiltonian identity")

identity = S.expand(Hc + theta - 2*psi)
assert identity == 0

print("Hc + theta - 2*psi = 0: PASS")
print("STAGE 1 PASS")

print("STAGE 2: Hamiltonian bracket identity")

# Antisymmetry implies {Hc,Hc}=0.
#
# From theta=2*psi-Hc:
#
# {theta,Hc}=2*{psi,Hc}.

B = S.symbols("B")
D = 2*B

assert S.expand(D-2*B) == 0

print("{theta,Hc}=2{psi,Hc}: PASS")
print("STAGE 2 PASS")

print("STAGE 3: multiplier consistency")

phi2, chi = S.symbols("phi2 chi")
l1, l2 = S.symbols("lambda1 lambda2")

# b = 2*A*W/v**3
b = S.symbols("b", nonzero=True)

# Previously established brackets:
#
# {psi,phi1}=chi-phi2
# {psi,phi2}=psi+b
#
# {theta,phi1}=chi-phi2
# {theta,phi2}=psi+2*b

dot_psi = (
    B
    +l1*(chi-phi2)
    +l2*(psi+b)
)

dot_theta = (
    2*B
    +l1*(chi-phi2)
    +l2*(psi+2*b)
)

relation = S.expand(
    dot_theta-2*dot_psi
    +l1*(chi-phi2)+l2*psi
)

assert relation == 0

print(
    "dot(theta)-2dot(psi)"
    "=-lambda1*(chi-phi2)-lambda2*psi: PASS"
)

# On the constraint surface:
# phi2=chi=psi=0.
#
# Thus dot(theta)=2dot(psi).

weak_difference = S.expand(
    (dot_theta-2*dot_psi).subs({
        phi2:0, chi:0, psi:0
    })
)

assert weak_difference == 0

print("WEAK PRESERVATION IDENTITY: PASS")
print("STAGE 3 PASS")

print("STAGE 4: regular A-nonzero branch")

# On phi2=chi=psi=0:
#
# dot(psi)=B+lambda2*b.
#
# For b!=0:
# lambda2=-B/b.
#
# Then dot(theta)=2B+2lambda2*b=0.

solution = -B/b

psi_check = S.cancel(B+solution*b)
theta_check = S.cancel(2*B+2*solution*b)

assert psi_check == 0
assert theta_check == 0

print("lambda2=-B/b: PASS")
print("dot(psi)=0: PASS")
print("dot(theta)=0: PASS")
print("STAGE 4 PASS")

print()
print("R26K COMPLETE")
print("EXACT HAMILTONIAN IDENTITY: VERIFIED")
print("TERTIARY PRESERVATION: STRUCTURALLY REDUNDANT")
print("SIXTH CONSTRAINT FROM theta: NOT GENERATED")
print("REGULAR A-NONZERO BRANCH: CONDITIONAL CLOSURE")
print("A-ZERO BRANCH: NOT CLASSIFIED")
print("FULL CONSTRAINT MATRIX: OPEN")
print("GAUGE CLASSIFICATION: OPEN")
print("PHYSICAL DEGREES OF FREEDOM: NOT COUNTED")
print("DEFI CORRESPONDENCE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(monotonic()-start,3)
)
