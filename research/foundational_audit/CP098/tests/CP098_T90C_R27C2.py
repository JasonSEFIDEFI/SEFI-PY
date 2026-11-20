import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27C2 ===", flush=True)

# Reduced scalar representation of the exact B expression.
# The scalar variables represent contractions in 4D.
#
# Primary constraint generators:
# phi1 = u.P
# phi2 = q.P
#
# Canonical transformations:
# {u,phi1}=0, {q,phi1}=u, {p1,phi1}=-P
# {u,phi2}=0, {q,phi2}=q, {P,phi2}=-P

s, v, d, W = S.symbols("s v d W", positive=True)
A, C, M, J = S.symbols("A C M J", real=True)
e, f, g, P2 = S.symbols("e f g P2", real=True)
phi1, phi2 = S.symbols("phi1 phi2", real=True)

# e = p0.q
# f = p1.eta.P
# g = p0.u
# P2 = P.eta.P
# J = sum(A_I pI/hI)

B = (
    -e
    -v*W*f/(2*C)
    +d*W*P2/(4*C*v)
    +5*A*d*W/v**5
    +d*M/v
    -2*A*W/(C*v**2)*(phi2+d*phi1/s)
    +2*W*J/v**2
)

print("STAGE 1: define generator actions", flush=True)

# Action of phi1 on scalar variables.
# {F,phi1} is implemented as a derivation.
#
# q -> u
# p1 -> -P
# d -> u.eta.u = -s
# W -> 0
# e -> g
# f -> -P2
# phi2 -> phi1

D1 = {
    s:0, v:0, d:-s, W:0,
    A:0, C:0, M:0, J:0,
    e:g, f:-P2, g:0, P2:0,
    phi1:0, phi2:phi1
}

# Action of phi2:
# q -> q, P -> -P
# d -> d, W -> 2W
# e -> e, f -> -f, P2 -> -2P2
# phi1 -> -phi1, phi2 -> 0

D2 = {
    s:0, v:0, d:d, W:2*W,
    A:0, C:0, M:0, J:0,
    e:e, f:-f, g:0, P2:-2*P2,
    phi1:-phi1, phi2:0
}

def derivation(expr, rules):
    return S.expand(
        sum(S.diff(expr,x)*dx for x,dx in rules.items())
    )

def check(name, expr):
    result = S.factor(expr)
    if result != 0:
        raise AssertionError(f"{name}: {result}")
    print("PASS:", name, flush=True)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: derive brackets", flush=True)

br1 = derivation(B,D1)
br2 = derivation(B,D2)

# Compare with compact candidate expressions.
expected1 = (
    -g
    +v*W*P2/(2*C)
    -s*W*P2/(4*C*v)
    -5*A*s*W/v**5
    -s*M/v
    -2*A*W/(C*v**2)*(
        phi1 - phi1 - d*phi1/s
    )
)

# We do not assume expected1 is correct.
# Print the derived expressions instead.

print("Bracket1 =", S.factor(br1), flush=True)
print("Bracket2 =", S.factor(br2), flush=True)

print("STAGE 2 PASS", flush=True)

print("STAGE 3: test R27A candidates", flush=True)

delta1 = S.factor(br1-B)
delta2 = S.factor(br2-2*B)

print("Delta1 =", delta1, flush=True)
print("Delta2 =", delta2, flush=True)

print("STAGE 3 PASS", flush=True)

print()
print("R27C2 COMPLETE")
print("SCALAR GENERATOR BRACKETS: DERIVED")
print("FULL PHASE-SPACE CROSS-CHECK: OPEN")
print("SINGULAR B=0 PRESERVATION: OPEN")
print("DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",round(perf_counter()-start,3))
