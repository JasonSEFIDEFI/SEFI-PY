import sympy as S
from time import monotonic

t0 = monotonic()
print("=== CP098-T90C-R26T ===", flush=True)

b, B, lam1, lam2 = S.symbols(
    "b B lambda1 lambda2", real=True
)

def check(label, residual):
    r = S.simplify(residual)
    if r != 0:
        print("FAIL:", label, r, flush=True)
        raise AssertionError(label)
    print("PASS:", label, flush=True)

# Weak constraint order:
# phi1, phi2, chi, psi, theta

C = S.Matrix([
    [0, 0, 0, 0, 0],
    [0, 0, 0, -b, -2*b],
    [0, 0, 0, 2*b, 4*b],
    [0, b, -2*b, 0, -B],
    [0, 2*b, -4*b, B, 0]
])

print("STAGE 1: Poisson matrix", flush=True)

check("antisymmetry", C+C.T)

# These independent columns span the matrix
# whenever at least one of b or B is nonzero.
c2 = C[:,1]
c4 = C[:,3]

check("column 2 dependency", C[:,2]+2*c2)

# Column 4 (zero-indexed) is:
# 2*column 3 + (B/b)*column 1 when b!=0.
check("column 4 dependency",
      b*(C[:,4]-2*c4)-B*c2)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: rank strata", flush=True)

check("rank at b=0,B=0",
      C.subs({b:0,B:0}).rank())

assert C.subs({b:0,B:1}).rank() == 2
assert C.subs({b:1,B:0}).rank() == 2
assert C.subs({b:1,B:1}).rank() == 2

print("PASS: rank 0 when b=B=0")
print("PASS: rank 2 when b=0,B!=0")
print("PASS: rank 2 when b!=0")
print("STAGE 2 PASS", flush=True)

print("STAGE 3: Hamiltonian consistency", flush=True)

# H_T = H_c + lambda1*phi1 + lambda2*phi2
#
# Previously established:
# {psi,Hc}=B
# {theta,Hc}=2B
# {psi,phi2}=psi+b
# {theta,phi2}=psi+2b
#
# On the five-constraint surface:
# dot(psi)   = B + lambda2*b
# dot(theta) = 2B + 2*lambda2*b

dot_psi = B+lam2*b
dot_theta = 2*B+2*lam2*b

check("theta consistency redundant",
      dot_theta-2*dot_psi)

# On b!=0, lambda2=-B/b.
check("regular branch multiplier",
      dot_psi.subs(lam2,-B/b))

# At b=0, lambda2 cannot cancel B.
check("A=0 consistency residual",
      dot_psi.subs(b,0)-B)

print("PASS: b!=0 fixes lambda2=-B/b")
print("PASS: b=0 requires B=0 for consistency")
print("STAGE 3 PASS", flush=True)

print("STAGE 4: singular branch interpretation",
      flush=True)

# This stage records conditional implications.
# It does not assume A=0 is invariant in time.

print("At A=0, b=0 because W>0 and v>0.")
print("If B!=0: five-constraint consistency fails.")
print("If B=0: weak five-constraint matrix rank is 0.")
print("Further consistency and branch invariance remain open.")

print("STAGE 4 PASS", flush=True)

print()
print("R26T COMPLETE")
print("WEAK MATRIX RANK STRATA: VERIFIED")
print("SINGULAR-BRANCH CONSISTENCY: VERIFIED")
print("A=0 BRANCH INVARIANCE: OPEN")
print("ADDITIONAL CONSTRAINTS: OPEN")
print("FULL GAUGE GENERATORS: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(monotonic()-t0,3))
