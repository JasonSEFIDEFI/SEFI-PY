"""CP098-T89A: constraint-bracket matrix construction.

Exact canonical brackets at a regular adapted point.
Five constraints, ten independent bracket entries.
No rank classification claimed.
"""
import sympy as sp

print("=== CP098-T89A ===", flush=True)

eta = sp.diag(-1,1,1,1)
eps = sp.symbols("eps", real=True)

def dot(x,y):
    return (x.T*eta*y)[0]

def pair(p,x):
    return (p.T*x)[0]

I,pI = sp.symbols("I pI",real=True)

A = sp.Rational(3,2)+I+I**2/3
C = sp.Rational(5,3)+I/2+I**2/5
B = (
    sp.Rational(7,4)
    +sp.Rational(2,5)*pI
    +I*pI/3+pI**2/7+I**2/11
)

def constraints(z):
    u = z[0:4,0]
    q = z[4:8,0]
    P = z[8:12,0]
    j = z[12:16,0]
    p0 = z[16:20,0]
    ii = z[20]
    pp = z[21]

    aa = A.subs(I,ii)
    cc = C.subs(I,ii)
    bb = B.subs({I:ii,pI:pp})

    v = sp.sqrt(-dot(u,u))
    T = u/v
    w = q+T*dot(T,q)
    L = sp.sqrt(dot(w,w))
    N = w/L
    k = L/v**2

    phiT = pair(P,T)
    phiN = pair(P,N)
    chiT = pair(j,T)
    chiN = pair(j,N)+2*aa*k/v
    psiT = (
        pair(p0,T)+bb-aa*k**2
        -v**4*k**2*dot(P,P)/(4*cc)
    )

    return (phiT,phiN,chiT,chiN,psiT)

# Canonical ordering:
# (u,j), (q,P), (I,pI).
# p0 is included as a fixed covector for these brackets.
J = sp.zeros(22)

for i in range(4):
    J[i,12+i] = 1
    J[12+i,i] = -1
    J[4+i,8+i] = 1
    J[8+i,4+i] = -1

J[20,21] = 1
J[21,20] = -1

assert J+J.T == sp.zeros(22)

# Adapted point:
# v=4, k=1/2, q time component=3.
# P and p1 satisfy the four primary/secondary
# constraints. p0 satisfies psiT.
v = sp.Integer(4)
k = sp.Rational(1,2)
aa = sp.Rational(3,2)
cc = sp.Rational(5,3)
bb = sp.Rational(7,4)

P = sp.Matrix([0,0,2,3])
j = sp.Matrix([0,-2*aa*k/v,5,7])
R = dot(P,P)

p0T = aa*k**2+v**4*k**2*R/(4*cc)-bb
p0 = sp.Matrix([p0T,11,13,17])

z0 = sp.Matrix(
    [4,0,0,0]
    +[3,8,0,0]
    +list(P)
    +list(j)
    +list(p0)
    +[0,0]
)

values = [
    sp.simplify(f)
    for f in constraints(z0)
]

for i,value in enumerate(values):
    print(f"CONSTRAINT_{i+1} =",value,flush=True)
    assert value == 0, (
        f"STOP: constraint {i+1} not satisfied"
    )

def gradient(index):
    result = []
    for i in range(22):
        direction = sp.eye(22)[:,i]
        expression = constraints(
            z0+eps*direction
        )[index]
        result.append(
            sp.diff(expression,eps).subs(eps,0)
        )
    return sp.Matrix(result)

print("Computing exact constraint gradients...",
      flush=True)

G = [gradient(i) for i in range(5)]

M = sp.zeros(5)

for i in range(5):
    for j in range(i+1,5):
        M[i,j] = sp.simplify(
            (G[i].T*J*G[j])[0]
        )
        M[j,i] = -M[i,j]

assert M+M.T == sp.zeros(5)

print("CONSTRAINT_BRACKET_MATRIX =",flush=True)
sp.pprint(M)

print("PASS: all ten independent entries evaluated")
print("PASS: constraint surface verified")
print("RANK CLASSIFICATION: NOT YET PERFORMED")
print("GENERAL-POINT MATRIX: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T89A COMPLETE")
