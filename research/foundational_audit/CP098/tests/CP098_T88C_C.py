"""CP098-T88C-C: exact nested Jacobi audit.

Triple (phi_N, chi_N, psi_T).
Full canonical phase space, explicit internal functions.
One non-adapted point; no general-point proof.
"""
import sympy as sp

print("=== CP098-T88C-C ===", flush=True)

eta = sp.diag(-1, 1, 1, 1)
e, f = sp.symbols("e f", real=True)

def dot(x, y):
    return (x.T*eta*y)[0]

def pair(p, x):
    return (p.T*x)[0]

def unpack(z):
    return (
        z[0:4, 0],
        z[4:8, 0],
        z[8:12, 0],
        z[12:16, 0],
        z[16],
        z[17]
    )

def constraints(z):
    u, q, P, j, I, pI = unpack(z)

    v = sp.sqrt(-dot(u,u))
    T = u/v
    w = q+T*dot(T,q)
    L = sp.sqrt(dot(w,w))
    N = w/L
    k = L/v**2

    A = sp.Rational(3,2)+I+I**2/3
    C = sp.Rational(5,3)+I/2+I**2/5
    B = (
        sp.Rational(7,4)
        +sp.Rational(2,5)*pI
        +I*pI/3+pI**2/7+I**2/11
    )

    phiN = pair(P,N)
    chiN = pair(j,N)+2*A*k/v

    # p0 is a fixed external covector for this triple.
    p0 = sp.Matrix([23,29,31,37])

    psiT = (
        pair(p0,T)+B-A*k**2
        -L**2*dot(P,P)/(4*C)
    )

    return phiN,chiN,psiT

z0 = sp.Matrix(
    [5,3,0,0]
    +[3,5,0,0]
    +[2,3,5,7]
    +[11,13,17,19]
    +[0,0]
)

# J maps a Hamiltonian gradient to its canonical
# Hamiltonian vector field.
J = sp.zeros(18)

for i in range(4):
    J[i,12+i] = 1
    J[12+i,i] = -1
    J[4+i,8+i] = 1
    J[8+i,4+i] = -1

J[16,17] = 1
J[17,16] = -1

assert J+J.T == sp.zeros(18)
assert J*J == -sp.eye(18)

names = ("phiN","chiN","psiT")

def directional_gradient(which):
    result = []
    for i in range(18):
        h = sp.eye(18)[:,i]
        expression = constraints(z0+e*h)[which]
        result.append(
            sp.simplify(
                sp.diff(expression,e).subs(e,0)
            )
        )
    return sp.Matrix(result)

print("Computing exact canonical gradients...",
      flush=True)

grad = [
    directional_gradient(i)
    for i in range(3)
]

X = [J*g for g in grad]

print("PASS: all canonical directions constructed",
      flush=True)

# For a scalar function F, compute its exact
# Hessian contraction with two fixed vectors.
# Both vectors are held constant while differentiating.
def hessian_contraction(which, a, b):
    expression = constraints(
        z0+e*a+f*b
    )[which]

    result = sp.diff(
        sp.diff(expression,e),f
    ).subs({e:0,f:0})

    return sp.simplify(result)

# {F,{G,H}} =
# -H_G(X_F,X_H)+H_H(X_F,X_G)
#
# All X vectors are evaluated at the reference
# point, but the underlying constraints are
# differentiated before evaluation.

cyclic = [(0,1,2),(1,2,0),(2,0,1)]
terms = []

for index,(i,j,k) in enumerate(cyclic,1):
    left = hessian_contraction(
        j,X[i],X[k]
    )
    right = hessian_contraction(
        k,X[i],X[j]
    )

    term = sp.simplify(-left+right)
    terms.append(term)

    print(
        f"CYCLIC_TERM_{index} =",
        term,flush=True
    )

residual = sp.simplify(sum(terms))

print("JACOBI_RESIDUAL =",residual,flush=True)
assert residual == 0, "STOP: Jacobi mismatch"

# Independent reverse-order second derivative
# check for the first pair of Hamiltonian directions.
forward = hessian_contraction(2,X[0],X[1])
reverse = hessian_contraction(2,X[1],X[0])

symmetry = sp.simplify(forward-reverse)

print("HESSIAN_SYMMETRY_RESIDUAL =",
      symmetry,flush=True)

assert symmetry == 0

print("PASS: full canonical nested Jacobi audit")
print("EXACT NON-ADAPTED POINT: VERIFIED")
print("ARBITRARY-POINT PROOF: OPEN")
print("GENERAL INTERNAL-FUNCTION PROOF: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T88C-C COMPLETE")
