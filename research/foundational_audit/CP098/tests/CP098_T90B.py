"""CP098-T90B: symbolic total-Hamiltonian tangency.

Adapted regular family, variable geometry, momenta,
and internal first-order data.

Does not assume tangency compatibility.
"""
import sympy as sp

print("=== CP098-T90B ===", flush=True)

Q = sp.Rational
eta = sp.diag(-1,1,1,1)

v0,k0 = sp.symbols(
    "v0 k0", positive=True, finite=True
)
a,r,s,x,y,n = sp.symbols(
    "a r s x y n", real=True
)
alpha,gamma = sp.symbols(
    "alpha gamma", real=True, nonzero=True
)
beta,ad,cd,bp = sp.symbols(
    "beta ad cd bp", real=True
)
lamT,lamN = sp.symbols(
    "lambda_T lambda_N", real=True
)

u = sp.Matrix(sp.symbols("u0:4", real=True))
q = sp.Matrix(sp.symbols("q0:4", real=True))
P = sp.Matrix(sp.symbols("P0:4", real=True))
j = sp.Matrix(sp.symbols("j0:4", real=True))
p0 = sp.Matrix(sp.symbols("z0:4", real=True))
I,pI = sp.symbols("I pI", real=True)

coords = list(u)+list(q)+list(P)+list(j)+list(p0)+[I,pI]

assert len(coords)==22

def dot(X,Y):
    return (X.T*eta*Y)[0]

def pair(X,Y):
    return (X.T*Y)[0]

A = alpha+ad*I
C = gamma+cd*I
B = beta+bp*pI

v = sp.sqrt(-dot(u,u))
T = u/v
w = q+T*dot(T,q)
L = sp.sqrt(dot(w,w))
k = L/v**2
N = w/L

F = [
    pair(P,T),
    pair(P,N),
    pair(j,T),
    pair(j,N)+2*A*k/v,
    pair(p0,T)+B-A*k**2
        -L**2*dot(P,P)/(4*C)
]

H = (
    pair(p0,u)+pair(j,q)
    +v*(B+A*k**2)
    -v**5*k**2*dot(P,P)/(4*C)
)

p0T = (
    alpha*k0**2
    +v0**4*k0**2*(r**2+s**2)/(4*gamma)
    -beta
)

surface = dict(zip(u,[v0,0,0,0]))
surface.update(dict(zip(q,[a,v0**2*k0,0,0])))
surface.update(dict(zip(P,[0,0,r,s])))
surface.update(dict(zip(j,[0,-2*alpha*k0/v0,x,y])))
surface.update(dict(zip(p0,[p0T,n,13,17])))
surface.update({I:0,pI:0})

print("Checking five constraints...",flush=True)

for index,f in enumerate(F,1):
    residual = sp.factor(f.subs(surface))
    print(f"CONSTRAINT_{index}_RESIDUAL =",
          residual,flush=True)
    assert residual == 0

# Canonical order:
# (u,j), (q,P), (I,pI).
J = sp.zeros(22)

for t in range(4):
    J[t,12+t]=1
    J[12+t,t]=-1
    J[4+t,8+t]=1
    J[8+t,4+t]=-1

J[20,21]=1
J[21,20]=-1

assert J+J.T == sp.zeros(22)

def gradient(expr,label):
    entries=[]
    for z in coords:
        entries.append(
            sp.factor(sp.diff(expr,z).subs(surface))
        )
    print(label+"_COMPLETE",flush=True)
    return sp.Matrix(entries)

G = [
    gradient(f,f"CONSTRAINT_GRADIENT_{i}")
    for i,f in enumerate(F,1)
]

GH = gradient(H,"HAMILTONIAN_GRADIENT")

print("Computing bracket matrix...",flush=True)

M=sp.zeros(5)

for i in range(5):
    for j2 in range(i+1,5):
        value=sp.factor((G[i].T*J*G[j2])[0])
        M[i,j2]=value
        M[j2,i]=-value

expected = {
    (1,3):-2*alpha/v0**3,
    (1,4):2*alpha*k0/v0**2,
    (2,3):4*alpha*k0/v0**2,
    (2,4):-4*alpha*k0**2/v0,
}

for (i,j2),value in expected.items():
    residual=sp.factor(M[i,j2]-value)
    print(f"M_{i+1}_{j2+1}_RESIDUAL =",
          residual,flush=True)
    assert residual==0

for j2 in range(5):
    assert sp.factor(M[0,j2])==0

assert sp.factor(M[1,2])==0

D=sp.factor(M[3,4])

print("D_BRACKET =",D,flush=True)

h=sp.Matrix([
    sp.factor((g.T*J*GH)[0])
    for g in G
])

for i,expr in enumerate(h,1):
    print(f"H_BRACKET_{i} =",expr,flush=True)

# The first three rows contain no primary
# multiplier. The last two rows must be
# mutually compatible.
obstructions=[
    sp.factor(h[0]),
    sp.factor(h[1]),
    sp.factor(h[2]),
    sp.factor(h[4]+k0*v0*h[3])
]

for i,expr in enumerate(obstructions,1):
    print(f"CONSISTENCY_OBSTRUCTION_{i} =",
          expr,flush=True)

# The fourth row determines lambda_N,
# provided alpha != 0 and v0 > 0.
lambdaN_solution=sp.factor(
    -h[3]/M[3,1]
)

print("LAMBDA_N_CANDIDATE =",
      lambdaN_solution,flush=True)

tangency=sp.Matrix([
    sp.factor(h[i]+lamT*M[i,0]+lamN*M[i,1])
    for i in range(5)
])

after=tangency.subs(
    lamN,lambdaN_solution
).applyfunc(sp.factor)

for i,expr in enumerate(after,1):
    print(f"POST_MULTIPLIER_RESIDUAL_{i} =",
          expr,flush=True)

assert after[3]==0

# Exact numerical reproduction of T90A.
case={
    v0:4,k0:Q(1,2),
    a:3,r:2,s:3,x:5,y:7,n:11,
    alpha:Q(3,2),gamma:Q(5,3),
    beta:Q(7,4),
    ad:1,cd:Q(1,2),bp:Q(2,5)
}

expected_h=sp.Matrix([
    0,0,0,-Q(64423,160),Q(64423,80)
])

numerical_residual=(
    h.subs(case)-expected_h
).applyfunc(sp.factor)

print("T90A_H_REPRODUCTION_RESIDUAL =",
      list(numerical_residual),flush=True)

assert numerical_residual==sp.zeros(5,1)

numerical_lambda=sp.factor(
    lambdaN_solution.subs(case)
)

print("T90A_LAMBDA_N =",numerical_lambda,
      flush=True)

assert numerical_lambda==Q(128846,15)

if all(expr==0 for expr in obstructions):
    print(
        "PASS: symbolic primary-multiplier "
        "compatibility on tested family",
        flush=True
    )
else:
    print(
        "STOP-GATE: nontrivial consistency "
        "conditions detected",
        flush=True
    )
    print(
        "DO NOT infer global Dirac closure",
        flush=True
    )

print("PASS: exact T90A reproduction")
print("GENERAL NON-ADAPTED CONSISTENCY: OPEN")
print("COMPLETE DIRAC ALGORITHM: OPEN")
print("GAUGE INTERPRETATION: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T90B COMPLETE")
