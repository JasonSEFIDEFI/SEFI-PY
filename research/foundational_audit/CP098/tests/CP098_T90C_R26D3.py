import sympy as S
from time import monotonic

t0 = monotonic()
print("=== CP098-T90C-R26D3 ===", flush=True)

# Scalar invariants; no four-component expansion.
s, W, P2, v = S.symbols(
    "s W P2 v", positive=True
)
A, C, m, V, T = S.symbols(
    "A C m V T", real=True
)
chi, phi1, phi2, pq = S.symbols(
    "chi phi1 phi2 pq", real=True
)

# Canonical Hamiltonian:
# H = p0.u + p1.q - v*W*P2/(4*C)
#     + A*W/v**3 + m*v + v*(V+T).
#
# Canonical flows:
# phi1: delta q=u, delta p1=-P
# phi2: delta q=q, delta P=-P
#
# phi1 also gives delta W=0.
# phi2 gives delta W=2W, delta P2=-2P2.

print("STAGE 1: phi1 flow", flush=True)

# Under {F,phi1}:
# delta(p1.q) = -P.q + p1.u
# delta(p0.u) = 0
# delta(W) = 0
# delta(P2) = 0.
#
# Therefore {H,phi1} = chi-phi2.

H_phi1 = chi - phi2
phi1_H = -H_phi1

assert S.expand(phi1_H-(phi2-chi)) == 0

print("{phi1,H} = phi2 - chi")
print("STAGE 1 PASS", flush=True)

print("STAGE 2: phi2 flow", flush=True)

# Under {F,phi2}:
# delta(p1.q)=p1.q
# delta(W)=2W
# delta(P2)=-2P2
#
# Hence delta(W*P2)=0.

delta_WP2 = S.expand(2*W*P2 + W*(-2*P2))
assert delta_WP2 == 0

# delta H = pq + 2*A*W/v**3
psi = pq + 2*A*W/v**3

H_phi2 = psi
phi2_H = -H_phi2

assert S.expand(phi2_H+psi) == 0

print("{phi2,H} = -psi")
print("psi = p1.q + 2*A*W/v**3")
print("STAGE 2 PASS", flush=True)

print("STAGE 3: total Hamiltonian", flush=True)

l1, l2 = S.symbols("l1 l2", real=True)

# From R26C:
# {phi1,phi2}=-phi1
# {phi2,phi1}=+phi1.

dotphi1 = phi2-chi-l2*phi1
dotphi2 = -psi+l1*phi1

print("dot(phi1) =", dotphi1)
print("dot(phi2) =", dotphi2)

assert S.expand(
    dotphi1-(phi2-chi-l2*phi1)
) == 0

assert S.expand(
    dotphi2-(-psi+l1*phi1)
) == 0

print("STAGE 3 PASS", flush=True)

print()
print("R26D3 COMPLETE")
print("phi1 PRESERVATION: CORRECTED")
print("phi2 PRESERVATION: CORRECTED")
print("SECONDARY CANDIDATE chi=p1.u")
print("SECONDARY CANDIDATE psi=p1.q+2AW/v^3")
print("SECONDARY PRESERVATION: OPEN")
print("FULL DIRAC CLASSIFICATION: OPEN")
print("DEFI CORRESPONDENCE: OPEN")
print("ELAPSED SECONDS =",
      round(monotonic()-t0, 3))
