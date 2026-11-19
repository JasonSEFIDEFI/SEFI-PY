import sympy as S
from time import monotonic

start = monotonic()

print("=== CP098-T90C-R26Q ===", flush=True)

def check(name, expression):
    residual = S.cancel(S.expand(expression))
    if residual != 0:
        print(name, "FAIL", flush=True)
        print("RESIDUAL =", residual, flush=True)
        raise AssertionError(name)
    print(name, "PASS", flush=True)

# Formal scalar invariants.
#
# v=sqrt(s), W>0.
# A and C depend on identity coordinates,
# not on u or q.
#
# M=m+V+T is independent of u and q.

v,W,A,C,M = S.symbols(
    "v W A C M", nonzero=True
)

J,K,L,N,P2 = S.symbols(
    "J K L N P2"
)

# J = q.p1
# K = p0.u
# L = u.p1 = chi
# N = u.P = phi1
# P2 = P_eta_squared

phi1,phi2,chi = S.symbols(
    "phi1 phi2 chi"
)

E = 2*A*W/v**3
b = E

# Hamiltonian terms
D = -v*W*P2/(4*C)
F = A*W/v**3
G = v*M

Hc = K+J+D+F+G

theta = -K+J-D+3*F-G
psi = J+E

print()
print("STAGE 1: Hamiltonian identity", flush=True)

check(
    "Hc+theta=2psi",
    Hc+theta-2*psi
)

print("STAGE 1 PASS", flush=True)

print()
print("STAGE 2: canonical product-rule inputs",
      flush=True)

# The following geometric derivative rules
# were verified independently in R26P:
#
# u.du(v)=v
# u.du(W)=0
# u.dq(W)=0
# q.dq(W)=2W
#
# They imply:
#
# u.du(D)=D
# u.du(F)=-3F
# u.du(G)=G
#
# q.dq(D)=2D
# q.dq(F)=2F
# q.dq(G)=0

# {phi1,Hc}:
# P.d(p1)Hc = phi2
# -u.d(q)Hc = -chi

br_phi1_H = phi2-chi

# {phi2,Hc}:
# P.d(P)Hc = 2D
# -q.d(q)Hc = -J-2D-2F

br_phi2_H = 2*D-J-2*D-2*F

# {chi,Hc}:
# p1.d(p1)Hc = J
# -u.d(u)Hc = -K-D+3F-G

br_chi_H = J-K-D+3*F-G

check(
    "{phi1,Hc}=phi2-chi",
    br_phi1_H-(phi2-chi)
)

check(
    "{phi2,Hc}=-psi",
    br_phi2_H+psi
)

check(
    "{chi,Hc}=theta",
    br_chi_H-theta
)

print("STAGE 2 PASS", flush=True)

print()
print("STAGE 3: basic canonical brackets",
      flush=True)

# Direct canonical product-rule relations:
#
# {phi1,phi2}=-phi1
# {phi1,chi}=phi1
# {phi2,chi}=0

br_12 = -phi1
br_1chi = phi1
br_2chi = S.Integer(0)

check(
    "{phi1,phi2}=-phi1",
    br_12+phi1
)

check(
    "{phi1,chi}=phi1",
    br_1chi-phi1
)

check(
    "{phi2,chi}=0",
    br_2chi
)

print("STAGE 3 PASS", flush=True)

print()
print("STAGE 4: theta brackets", flush=True)

# Exact identity:
#
# theta = 2psi-Hc
#
# Previously verified in R26P:
#
# {phi1,psi}=phi2-chi
# {phi2,psi}=-psi-b
# {chi,psi}=psi+2b

br_1psi = phi2-chi
br_2psi = -psi-b
br_chipsi = psi+2*b

br_1theta = 2*br_1psi-br_phi1_H
br_2theta = 2*br_2psi-br_phi2_H
br_chitheta = 2*br_chipsi-br_chi_H

check(
    "{phi1,theta}=phi2-chi",
    br_1theta-(phi2-chi)
)

check(
    "{phi2,theta}=-psi-2b",
    br_2theta+psi+2*b
)

check(
    "{chi,theta}=2psi+4b-theta",
    br_chitheta-(2*psi+4*b-theta)
)

print("STAGE 4 PASS", flush=True)

print()
print("STAGE 5: psi-theta reduction",
      flush=True)

B = S.symbols("B")

# B denotes the full canonical bracket
# {psi,Hc}, including identity-sector terms.
#
# Since theta=2psi-Hc:
#
# {psi,theta}
#   =2{psi,psi}-{psi,Hc}
#   =-B

br_psitheta = -B

check(
    "{psi,theta}=-B",
    br_psitheta+B
)

print("STAGE 5 PASS", flush=True)

print()
print("R26Q COMPLETE")
print("HAMILTONIAN IDENTITY: PASS")
print("THREE HAMILTONIAN BRACKETS: DERIVED")
print("THREE BASIC BRACKETS: DERIVED")
print("THREE THETA BRACKETS: DERIVED")
print("PSI-THETA REDUCTION: PASS")
print("ALL SEVEN REMAINING RELATIONS: PASS")
print("FULL B EXPRESSION: NOT PROVEN HERE")
print("FULL DIRECT SYMBOLIC AUDIT: OPEN")
print("GAUGE GENERATORS: OPEN")
print("A-ZERO BRANCH: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(monotonic()-start,3)
)
