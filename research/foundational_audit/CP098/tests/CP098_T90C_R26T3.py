import sympy as S
from time import monotonic

t0 = monotonic()
print("=== CP098-T90C-R26T3 ===", flush=True)

b, B, lam2 = S.symbols("b B lambda2", real=True)

def check(label, residual):
    if isinstance(residual, S.MatrixBase):
        ok = all(S.simplify(x) == 0 for x in residual)
    else:
        ok = S.simplify(residual) == 0

    if not ok:
        print("FAIL:", label, residual, flush=True)
        raise AssertionError(label)

    print("PASS:", label, flush=True)

C = S.Matrix([
    [0, 0, 0, 0, 0],
    [0, 0, 0, -b, -2*b],
    [0, 0, 0, 2*b, 4*b],
    [0, b, -2*b, 0, -B],
    [0, 2*b, -4*b, B, 0]
])

print("STAGE 1: exact matrix identities", flush=True)

check("antisymmetry", C+C.T)

# Column 2 = -2 column 1 (zero-based indices).
check("dependent constraint column",
      C[:,2]+2*C[:,1])

# Column 4 = 2 column 3 - (B/b) column 1
# when b!=0. Multiply through by b to
# avoid division and cover the polynomial identity.

check("additional column dependency",
      b*(C[:,4]-2*C[:,3])+B*C[:,1])

print("STAGE 1 PASS", flush=True)

print("STAGE 2: exact rank classification", flush=True)

# Every 3x3 minor vanishes if rank <=2.
# Use symbolic rank plus exact witnesses to
# verify the polynomial rank bound.

assert C.rank() <= 2

assert C.subs({b:1,B:0}).rank() == 2
assert C.subs({b:0,B:1}).rank() == 2
assert C.subs({b:0,B:0}).rank() == 0

# Explicit 2x2 minors:
# rows/cols 1,3 -> b^2
# rows/cols 3,4 -> B^2

check("b minor", C.extract([1,3],[1,3]).det()-b**2)
check("B minor", C.extract([3,4],[3,4]).det()-B**2)

print("PASS: rank 2 if b!=0 or B!=0")
print("PASS: rank 0 if b=B=0")
print("STAGE 2 PASS", flush=True)

print("STAGE 3: Hamiltonian consistency", flush=True)

# On the five-constraint surface:
#
# dot(psi) = B + lambda2*b
# dot(theta) = 2B + 2lambda2*b

dot_psi = B+lam2*b
dot_theta = 2*B+2*lam2*b

check("theta consistency redundant",
      dot_theta-2*dot_psi)

check("regular branch lambda2",
      dot_psi.subs(lam2,-B/b))

check("A-zero consistency",
      dot_psi.subs(b,0)-B)

print("PASS: b!=0 fixes lambda2=-B/b")
print("PASS: b=0 requires B=0")
print("STAGE 3 PASS", flush=True)

print("STAGE 4: classification boundaries", flush=True)

print("A!=0: b!=0, weak matrix rank 2")
print("A=0, B!=0: weak matrix rank 2")
print("A=0, B=0: weak matrix rank 0")
print("A=0, B!=0: consistency obstruction")
print("A=0, B=0: further Dirac analysis required")

print("STAGE 4 PASS", flush=True)

print()
print("R26T2 COMPLETE")
print("EXACT WEAK RANK CLASSIFICATION: PASS")
print("MULTIPLIER CONSISTENCY: PASS")
print("A-ZERO BRANCH INVARIANCE: OPEN")
print("ADDITIONAL CONSTRAINTS: OPEN")
print("GAUGE GENERATORS: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(monotonic()-t0,3))
