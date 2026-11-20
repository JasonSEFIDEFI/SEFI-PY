"""CP098-T89D: symbolic regular-family rank audit.

Adapted geometry with v=4, kappa=1/2.
Symbolic a,r,s,x,y,n.
Exact canonical differentiation before surface substitution.
"""
import sympy as sp

print("=== CP098-T89D ===", flush=True)

eta = sp.diag(-1,1,1,1)

a,r,s,x,y,n = sp.symbols(
    "a r s x y n", real=True
)

I,pI = sp.symbols("I pI", real=True)

A = sp.Rational(3,2)+I+I**2/3
C = sp.Rational(5,3)+I/2+I**2/5
B = (
    sp.Rational(7,4)
    +sp.Rational(2,5)*pI
    +I*pI/3+pI**2/7+I**2/11
)

def dot(X,Y):
    return (X.T*eta*Y)[0]

def pair(p,X):
    return (p.T*X)[0]

# Independent canonical coordinates.
u = sp.Matrix(sp.symbols("u0:4", real=True))
q = sp.Matrix(sp.symbols("q0:4", real=True))
P = sp.Matrix(sp.symbols("P0:4", real=True))
j = sp.Matrix(sp.symbols("j0:4", real=True))
p0 = sp.Matrix(sp.symbols("z0:4", real=True))

v = sp.sqrt(-dot(u,u))
T = u/v
w = q+T*dot(T,q)
L = sp.sqrt(dot(w,w))
N = w/L
k = L/v**2

phiT = pair(P,T)
phiN = pair(P,N)
chiT = pair(j,T)
chiN = pair(j,N)+2*A*k/v
psiT = (
    pair(p0,T)+B-A*k**2
    -L**2*dot(P,P)/(4*C)
)

F = [phiT,phiN,chiT,chiN,psiT]

# Canonical pairs (u,j), (q,P), (I,pI).
# The (y,p0) contribution is zero because
# none of these constraints depends on y.
coords = list(u)+list(q)+list(P)+list(j)+[I,pI]

J = sp.zeros(18)

for t in range(4):
    J[t,12+t] = 1
    J[12+t,t] = -1
    J[4+t,8+t] = 1
    J[8+t,4+t] = -1

J[16,17] = 1
J[17,16] = -1

assert J+J.T == sp.zeros(18)

v0 = sp.Integer(4)
k0 = sp.Rational(1,2)
A0 = sp.Rational(3,2)
C0 = sp.Rational(5,3)
B0 = sp.Rational(7,4)

R = r**2+s**2

p0T = (
    A0*k0**2
    +v0**4*k0**2*R/(4*C0)
    -B0
)

surface = {
    **dict(zip(u,[4,0,0,0])),
    **dict(zip(q,[a,8,0,0])),
    **dict(zip(P,[0,0,r,s])),
    **dict(zip(j,[0,-2*A0*k0/v0,x,y])),
    **dict(zip(p0,[p0T,n,13,17])),
    I:sp.S.Zero,
    pI:sp.S.Zero,
}

print("Verifying symbolic constraint surface...",
      flush=True)

for idx,f in enumerate(F,1):
    residual = sp.simplify(f.subs(surface))
    print(f"CONSTRAINT_{idx}_RESIDUAL =",
          residual,flush=True)
    assert residual == 0

# Compute derivatives before imposing constraints.
# This distinction is essential: differentiating
# a constraint already set to zero is invalid.
print("Computing symbolic canonical gradients...",
      flush=True)

gradients = []

for idx,f in enumerate(F,1):
    g = sp.Matrix([
        sp.diff(f,z).subs(surface)
        for z in coords
    ])
    gradients.append(g)
    print(f"GRADIENT_{idx}_COMPLETE",flush=True)

M = sp.zeros(5)

for i in range(5):
    for jdx in range(i+1,5):
        M[i,jdx] = sp.factor(
            (gradients[i].T*J*gradients[jdx])[0]
        )
        M[jdx,i] = -M[i,jdx]

print("SYMBOLIC_MATRIX =",flush=True)
sp.pprint(M)

assert M+M.T == sp.zeros(5)

K = M[1:5,1:5]

pf = sp.factor(
    K[0,1]*K[2,3]
    -K[0,2]*K[1,3]
    +K[0,3]*K[1,2]
)

print("SYMBOLIC_PFAFFIAN =",pf,flush=True)

# Independent nonzero minor.
minor = sp.factor(
    M.extract([1,3],[1,3]).det()
)

print("SYMBOLIC_2X2_MINOR =",minor,flush=True)

# Compare the symbolic family with both
# previously verified numerical matrices.
case1 = {a:3,r:2,s:3,x:5,y:7,n:11}

expected1 = sp.Matrix([
    [0,0,0,0,0],
    [0,0,0,-sp.Rational(3,64),sp.Rational(3,32)],
    [0,0,0,sp.Rational(3,16),-sp.Rational(3,8)],
    [0,sp.Rational(3,64),-sp.Rational(3,16),
     0,-sp.Rational(64333,640)],
    [0,-sp.Rational(3,32),sp.Rational(3,8),
     sp.Rational(64333,640),0],
])

comparison = sp.simplify(
    M.subs(case1)-expected1
)

print("T89A_REPRODUCTION_RESIDUAL =",
      comparison,flush=True)

assert comparison == sp.zeros(5)

if pf != 0:
    raise AssertionError(
        "STOP: symbolic Pfaffian is nonzero"
    )

if minor == 0:
    raise AssertionError(
        "STOP: rank-two lower bound missing"
    )

# First row and column must vanish
# for the rank-two inference from K.
assert all(
    sp.simplify(M[0,i]) == 0
    for i in range(5)
)

print("PASS: symbolic family has rank exactly two")
print("SYMBOLIC ADAPTED FAMILY: VERIFIED")
print("VARIABLE v AND kappa: OPEN")
print("GENERAL NON-ADAPTED SECTOR: OPEN")
print("FIRST-CLASS CLASSIFICATION: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T89D COMPLETE")
