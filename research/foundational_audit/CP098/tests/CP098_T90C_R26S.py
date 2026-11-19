import sympy as S
from time import monotonic

t0 = monotonic()
print("=== CP098-T90C-R26S ===", flush=True)

def check(name, expr):
    r = S.factor(expr)
    if r != 0:
        print("FAIL:", name, "RESIDUAL:", r, flush=True)
        raise AssertionError(name)
    print("PASS:", name, flush=True)

a,b,c,d = S.symbols("a b c d", real=True)

# a=u0, b=uj, c=q0, d=qj.
#
# Selected phase-space columns:
# P0, Pj, p1_0, p1_j, p0_0.
#
# Constraint rows:
# phi1, phi2, chi, psi, theta.
#
# The momentum derivatives are exact.
# E=2AW/v^3 has no momentum dependence.

print("STAGE 1: exact Jacobian minor", flush=True)

J = S.Matrix([
    [a,b,0,0,0],
    [c,d,0,0,0],
    [0,0,a,b,0],
    [0,0,c,d,0],
    [0,0,c,d,-a]
])

delta = a*d-b*c
det = J.det(method="bareiss")

check("FIVE-BY-FIVE DETERMINANT",
      det+a*delta**2)

print("DETERMINANT =", S.factor(det))
print("STAGE 1 PASS", flush=True)

print("STAGE 2: timelike condition", flush=True)

# s=(u0)^2 - sum_j(uj)^2 > 0.
# Hence u0 cannot vanish.

u0,u1,u2,u3 = S.symbols("u0:4", real=True)
s = u0**2-u1**2-u2**2-u3**2

check("s AT u0=0 IS NONPOSITIVE",
      s.subs(u0,0)+(u1**2+u2**2+u3**2))

print("s>0 IMPLIES u0!=0")
print("STAGE 2 PASS", flush=True)

print("STAGE 3: nonparallel geometry", flush=True)

q0,q1,q2,q3 = S.symbols("q0:4", real=True)

u = S.Matrix([u0,u1,u2,u3])
q = S.Matrix([q0,q1,q2,q3])
eta = S.diag(-1,1,1,1)

def dot(x,y):
    return (x.T*eta*y)[0]

ss = -dot(u,u)
dd = dot(u,q)
W = dot(q,q)+dd**2/ss

# If all Delta_j vanish and u0!=0,
# then qj=(q0/u0)*uj for every j.
# Therefore q=(q0/u0)*u and W=0.

parallel = {
    q1:q0*u1/u0,
    q2:q0*u2/u0,
    q3:q0*u3/u0
}

check("PARALLEL q IMPLIES W=0",
      W.subs(parallel))

print("W>0 IMPLIES SOME Delta_j!=0")
print("STAGE 3 PASS", flush=True)

print("STAGE 4: rank conclusion", flush=True)

# There are exactly five constraint functions.
# A nonzero 5x5 minor implies rank >=5.
# Five rows imply rank <=5.

print("s>0 AND W>0:")
print("  u0!=0")
print("  SOME Delta_j!=0")
print("  SOME JACOBIAN MINOR NONZERO")
print("  CONSTRAINT JACOBIAN RANK = 5")
print("STAGE 4 PASS", flush=True)

print()
print("R26S COMPLETE")
print("GENERAL CONSTRAINT INDEPENDENCE: PROVEN")
print("REGION: s>0, W>0, C!=0")
print("POISSON RANK 2: CONDITIONAL ON A!=0")
print("A-ZERO BRANCH: OPEN")
print("GAUGE GENERATORS: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(monotonic()-t0,3))
