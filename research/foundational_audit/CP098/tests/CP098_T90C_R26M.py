import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R26M ===", flush=True)

f1, f2, chi, psi, theta = S.symbols(
    "phi1 phi2 chi psi theta"
)
b, B = S.symbols("b B", real=True)

names = ["phi1", "phi2", "chi", "psi", "theta"]
constraints = [f1, f2, chi, psi, theta]

print("STAGE 1: canonical bracket inventory", flush=True)

# Previously derived identities.
# b = 2*A*W/v**3
# B = {psi,Hc}
#
# Upper triangular entries only.
# All theta entries are reconstructed
# using theta = 2*psi - Hc.

brackets = {
    (0,1): -f1,
    (0,2): f1,
    (0,3): f2-chi,
    (1,2): S.Integer(0),
    (1,3): -psi-b,
    (2,3): psi+2*b
}

# Known Hamiltonian brackets:
#
# {phi1,Hc}=phi2-chi
# {phi2,Hc}=-psi
# {chi,Hc}=theta
# {psi,Hc}=B

with_H = [
    f2-chi,
    -psi,
    theta,
    B
]

for i in range(4):
    brackets[(i,4)] = S.expand(
        2*brackets[(i,3)]-with_H[i]
    )

assert len(brackets)==10

for (i,j),value in sorted(brackets.items()):
    print(
        "{%s,%s} =" % (names[i],names[j]),
        value,
        flush=True
    )

print("STAGE 1 PASS", flush=True)

print("STAGE 2: construct matrix independently",
      flush=True)

C = S.zeros(5)

for (i,j),value in brackets.items():
    C[i,j] = value
    C[j,i] = -value

assert C+C.T == S.zeros(5)

weak_subs = {x:0 for x in constraints}
Cweak = C.subs(weak_subs)

expected = S.Matrix([
    [0, 0, 0, 0, 0],
    [0, 0, 0, -b, -2*b],
    [0, 0, 0, 2*b, 4*b],
    [0, b, -2*b, 0, -B],
    [0, 2*b, -4*b, B, 0]
])

assert Cweak == expected

print("ALL TEN WEAK ENTRIES: MATCH R26L2")
print("STAGE 2 PASS", flush=True)

print("STAGE 3: symbolic rank and kernel",
      flush=True)

assert Cweak.rank()==2

kernel = Cweak.nullspace()
assert len(kernel)==3

for vec in kernel:
    assert Cweak*vec == S.zeros(5,1)

print("RANK = 2")
print("NULLITY = 3")
print("STAGE 3 PASS", flush=True)

print("STAGE 4: explicit first-class combinations",
      flush=True)

# These are weakly first-class combinations
# of the five constraints.
#
# F1 = phi1
# F2 = 2*phi2+chi
# F3 = (B/b)*phi2-2*psi+theta
#
# Coefficients may depend on phase space.
# Derivatives of coefficients multiply
# constraints and therefore vanish weakly.

v1 = S.Matrix([1,0,0,0,0])
v2 = S.Matrix([0,2,1,0,0])
v3 = S.Matrix([0,B/b,0,-2,1])

for i,vec in enumerate([v1,v2,v3],1):
    residual = (Cweak*vec).applyfunc(S.cancel)
    assert residual == S.zeros(5,1)
    print("FIRST-CLASS COMBINATION",i,"PASS")

assert S.Matrix.hstack(v1,v2,v3).rank()==3

print("STAGE 4 PASS", flush=True)

print("STAGE 5: exact rational witness",
      flush=True)

witness = {
    b:S.Rational(387,512),
    B:S.Rational(1204251,81920)
}

Cw = Cweak.subs(witness)

assert Cw.rank()==2
assert len(Cw.nullspace())==3

print("WITNESS RANK =",Cw.rank())
print("WITNESS NULLITY =",len(Cw.nullspace()))
print("STAGE 5 PASS", flush=True)

print()
print("R26M COMPLETE")
print("TEN BRACKETS: ALGEBRAICALLY RECONSTRUCTED")
print("THETA BRACKETS: EXACT IDENTITY USED")
print("R26L2 WEAK MATRIX: REPRODUCED")
print("WEAK RANK: 2")
print("WEAK NULLITY: 3")
print("FIRST-CLASS COMBINATIONS: VERIFIED")
print("CANONICAL BRACKET PREMISES: PRIOR RESULTS")
print("GENERAL psi-Hc BRACKET: STILL OPEN")
print("A-NONZERO BRANCH: CONDITIONAL CLASSIFICATION")
print("A-ZERO BRANCH: OPEN")
print("GAUGE GENERATORS: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(monotonic()-start,3))
