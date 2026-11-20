import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R25B ===", flush=True)

lam = S.symbols("lam", real=True)
x = S.symbols("x", real=True)

# Explicit candidate identity trajectories.
IA = lam
IS = S.Integer(1)
IW = S.Integer(1)

# Explicit candidate coupling functions.
mu1 = 3 + 2*x
mu2 = 2 + x + S.Rational(3, 2)*x**2

A = S.expand(mu1.subs(x, IA)*IS)
C = S.expand(mu2.subs(x, IA)*IW)

def jet(f, order):
    return S.simplify(
        S.diff(f, lam, order).subs(lam, 0)
    )

print("STAGE 1: coupling jets", flush=True)

actual = (
    jet(A, 0), jet(A, 1),
    jet(C, 0), jet(C, 1), jet(C, 2)
)
expected = tuple(map(S.Integer, (3, 2, 2, 1, 3)))

print("ACTUAL =", actual)
print("EXPECTED =", expected)

assert actual == expected, "STOP: coupling-jet mismatch"
print("STAGE 1 PASS", flush=True)

print("STAGE 2: explicit chain rules", flush=True)

mu1_I = S.diff(mu1, x).subs(x, IA)
mu2_I = S.diff(mu2, x).subs(x, IA)
mu2_II = S.diff(mu2, x, 2).subs(x, IA)

Aprime = (
    mu1_I*S.diff(IA, lam)*IS
    + mu1.subs(x, IA)*S.diff(IS, lam)
)

Cprime = (
    mu2_I*S.diff(IA, lam)*IW
    + mu2.subs(x, IA)*S.diff(IW, lam)
)

Csecond = (
    mu2_II*S.diff(IA, lam)**2*IW
    + mu2_I*S.diff(IA, lam, 2)*IW
    + 2*mu2_I*S.diff(IA, lam)*S.diff(IW, lam)
    + mu2.subs(x, IA)*S.diff(IW, lam, 2)
)

assert S.simplify(S.diff(A, lam)-Aprime) == 0
assert S.simplify(S.diff(C, lam)-Cprime) == 0
assert S.simplify(S.diff(C, lam, 2)-Csecond) == 0

print("CHAIN-RULE RESIDUALS = [0, 0, 0]")
print("STAGE 2 PASS", flush=True)

print("STAGE 3: identity-source derivatives", flush=True)

# Geometric Lagrangian as a function of A and C.
ss, W, Z = S.symbols("ss W Z", positive=True)
Lgeom = (
    -A*W/ss**S.Rational(3, 2)
    -C*Z/(S.sqrt(ss)*W)
)

# Differentiation with respect to independent identity
# coordinates, before inserting the trajectory.
ia, is_, iw = S.symbols("ia is iw", real=True)
Afield = mu1.subs(x, ia)*is_
Cfield = mu2.subs(x, ia)*iw

Lfield = (
    -Afield*W/ss**S.Rational(3, 2)
    -Cfield*Z/(S.sqrt(ss)*W)
)

sources = [
    S.diff(Lfield, variable)
    for variable in (ia, is_, iw)
]

print("dLgeom/dIA =", sources[0])
print("dLgeom/dIS =", sources[1])
print("dLgeom/dIW =", sources[2])

# No identity-velocity dependence in Lgeom.
# Thus geometric identity canonical momenta vanish.
print("GEOMETRIC IDENTITY-VELOCITY MOMENTA = [0, 0, 0]")
print("STAGE 3 PASS", flush=True)

print()
print("R25B COMPLETE")
print("COUPLING-JET EXISTENCE: VERIFIED")
print("CHAIN-RULE CONSISTENCY: VERIFIED")
print("IDENTITY-SOURCE TERMS: DERIVED")
print("GENERAL COUPLING FUNCTIONS: NOT DETERMINED")
print("IDENTITY EQUATIONS OF MOTION: NOT SOLVED")
print("FULL HAMILTONIAN: OPEN")
print("DIRAC SYSTEM: OPEN")
print("DEFI CORRESPONDENCE: OPEN")
print("ELAPSED SECONDS =", round(monotonic()-start, 3))
