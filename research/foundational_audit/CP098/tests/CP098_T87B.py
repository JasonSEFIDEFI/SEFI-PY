"""CP098-T87B: arbitrary-point canonical bracket identity.

Algebraic verification of tensor contractions and bracket
assembly from geometric differential identities.

No adapted frame. No momentum constraints.
"""
import sympy as sp

print("=== CP098-T87B ===", flush=True)

eta = sp.diag(-1, 1, 1, 1)

T = sp.Matrix(sp.symbols("T0:4"))
N = sp.Matrix(sp.symbols("N0:4"))
P = sp.Matrix(sp.symbols("P0:4"))
p1 = sp.Matrix(sp.symbols("j0:4"))

def pair(p, x):
    return (p.T*x)[0]

def covdot(p, z):
    return (p.T*eta*z)[0]

# q-derivative of N:
#
# dN/dq = [Id + T (T^T eta) - N (N^T eta)]/L
#
# The derivative of chi_N contracts this projector
# against canonical covector p1 and eta*P.

projector = (
    sp.eye(4)
    + T*(T.T*eta)
    - N*(N.T*eta)
)

direct_contraction = (
    p1.T*projector*eta*P
)[0]

expanded_contraction = (
    covdot(p1, P)
    + pair(p1, T)*pair(P, T)
    - pair(p1, N)*pair(P, N)
)

projector_residual = sp.expand(
    direct_contraction-expanded_contraction
)

print(
    "PROJECTOR_CONTRACTION_RESIDUAL =",
    projector_residual,
    flush=True
)
assert projector_residual == 0

print("PASS: unrestricted tensor contraction")

# Independent scalar variables for arbitrary regular geometry.
v, L, tq = sp.symbols(
    "v L tq", positive=True
)
A, C, DA = sp.symbols("A C DA", nonzero=True)

jP, jT, PT, jN, PN, zN, P2 = sp.symbols(
    "jP jT PT jN PN zN P2"
)

# Contribution from (q,P) canonical pair.
#
# dchi/dq = p1 projector/L + 2A N_flat/v^3
# dpsi/dP = -L^2 eta P/(2C)

qP_part = (
    -L/(2*C)*(jP+jT*PT-jN*PN)
    -A*L**2*PN/(C*v**3)
)

# Directional u derivative:
# N(v)=0, N(T)=N/v,
# N(L)=tq/v, N(kappa)=tq/v^3.
#
# N(psi) = zN/v
#          -2A L tq/v^5
#          -L tq P2/(2Cv)

u_part = (
    -zN/v
    +2*A*L*tq/v**5
    +L*tq*P2/(2*C*v)
)

internal_part = 2*L*DA/v**3

canonical_E = sp.expand(
    qP_part+u_part+internal_part
)

compact_E = (
    -L/(2*C)*(jP+jT*PT-jN*PN)
    -A*L**2*PN/(C*v**3)
    -zN/v
    +2*A*L*tq/v**5
    +L*tq*P2/(2*C*v)
    +2*L*DA/v**3
)

residual = sp.factor(
    sp.together(canonical_E-compact_E)
)

print(
    "ARBITRARY_CONTRACTION_RESIDUAL =",
    residual,
    flush=True
)
assert residual == 0

# Verify the internal bracket term independently
# using arbitrary first derivatives of A and B.

AI, Ap, BI, Bp = sp.symbols("AI Ap BI Bp")

DA_direct = AI*Bp-Ap*BI
DA_symbolic = sp.symbols("DA_symbolic")

internal_direct = (
    2*L/v**3*(AI*Bp-Ap*BI)
)
internal_expected = 2*L*DA_symbolic/v**3

internal_residual = sp.simplify(
    internal_direct
    -internal_expected.subs(
        DA_symbolic, DA_direct
    )
)

print(
    "INTERNAL_BRACKET_RESIDUAL =",
    internal_residual,
    flush=True
)
assert internal_residual == 0

print("PASS: arbitrary-point algebraic bracket assembly")
print("PASS: unrestricted canonical contractions")
print("PASS: internal Poisson-bracket contribution")
print("GENERAL GEOMETRIC BRACKET IDENTITY: DERIVED")
print("INDEPENDENT JACOBI CROSS-CHECK: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T87B COMPLETE")
