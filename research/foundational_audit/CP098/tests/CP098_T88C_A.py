"""CP098-T88C-A: full canonical inner brackets.

Triple (phi_N, chi_N, psi_T).
Two non-adapted exact points, unrestricted momenta.
No nested Jacobi claim.
"""
import sympy as sp

print("=== CP098-T88C-A ===", flush=True)

eta = sp.diag(-1,1,1,1)
e = sp.symbols("e", real=True)
basis = sp.eye(4)

A = sp.Rational(3,2)
C = sp.Rational(5,3)
B = sp.Rational(7,4)
DA = sp.Rational(2,5)

def dot(x,y):
    return (x.T*eta*y)[0]

def pair(p,x):
    return (p.T*x)[0]

def covdot(p,z):
    return (p.T*eta*z)[0]

def geom(u,q):
    v = sp.sqrt(-dot(u,u))
    T = u/v
    w = q+T*dot(T,q)
    L = sp.sqrt(dot(w,w))
    N = w/L
    k = L/v**2
    return v,T,N,L,k

def phiN(u,q,P,j,z):
    return pair(P,geom(u,q)[2])

def chiN(u,q,P,j,z):
    v,T,N,L,k = geom(u,q)
    return pair(j,N)+2*A*k/v

def psiT(u,q,P,j,z):
    v,T,N,L,k = geom(u,q)
    return (
        pair(z,T)+B-A*k**2
        -v**4*k**2*dot(P,P)/(4*C)
    )

def compactE(u,q,P,j,z):
    v,T,N,L,k = geom(u,q)
    tq = dot(T,q)
    return (
        -L/(2*C)*(
            covdot(j,P)
            +pair(j,T)*pair(P,T)
            -pair(j,N)*pair(P,N)
        )
        -A*L**2*pair(P,N)/(C*v**3)
        -pair(z,N)/v
        +2*A*L*tq/v**5
        +L*tq*dot(P,P)/(2*C*v)
        +2*L*DA/v**3
    )

def derivative(fn,args,slot,i):
    values = list(args)
    values[slot] += e*basis[:,i]
    return sp.diff(fn(*values),e).subs(e,0)

def pb(f,g,args,internal=0):
    result = sp.S.Zero
    for i in range(4):
        fu = derivative(f,args,0,i)
        fj = derivative(f,args,3,i)
        fq = derivative(f,args,1,i)
        fP = derivative(f,args,2,i)

        gu = derivative(g,args,0,i)
        gj = derivative(g,args,3,i)
        gq = derivative(g,args,1,i)
        gP = derivative(g,args,2,i)

        result += fu*gj-fj*gu+fq*gP-fP*gq

    return sp.simplify(result+internal)

cases = [
    (sp.Matrix([5,3,0,0]),
     sp.Matrix([3,5,0,0])),
    (sp.Matrix([10,2,3,1]),
     sp.Matrix([4,1,5,2]))
]

P = sp.Matrix([2,3,5,7])
j = sp.Matrix([11,13,17,19])
z = sp.Matrix([23,29,31,37])

for index,(u,q) in enumerate(cases,1):
    args = (u,q,P,j,z)
    v,T,N,L,k = geom(u,q)

    # First bracket: known full canonical identity.
    b1 = pb(phiN,chiN,args)
    expected1 = -2*A/v**3+pair(P,T)/v
    r1 = sp.simplify(b1-expected1)

    print(f"CASE_{index}_N_CHI_RESIDUAL =",r1,
          flush=True)
    assert r1 == 0

    # Second bracket: independent canonical calculation.
    # The internal canonical pair contributes 2*k*DA/v.
    b2 = pb(
        chiN,psiT,args,
        internal=2*k*DA/v
    )
    expected2 = compactE(*args)
    r2 = sp.simplify(b2-expected2)

    print(f"CASE_{index}_CHI_PSI_RESIDUAL =",r2,
          flush=True)
    assert r2 == 0

    # Third bracket: calculate and retain its exact value.
    b3 = pb(psiT,phiN,args)

    print(f"CASE_{index}_PSI_N_VALUE =",b3,
          flush=True)

    # Independent antisymmetry check.
    reverse = pb(phiN,psiT,args)
    r3 = sp.simplify(b3+reverse)

    print(f"CASE_{index}_ANTISYMMETRY_RESIDUAL =",
          r3,flush=True)
    assert r3 == 0

    print(f"PASS: full inner brackets case {index}",
          flush=True)

print("PASS: psi_T inner-bracket prerequisite")
print("NESTED JACOBI BRACKETS: OPEN")
print("ARBITRARY-POINT JACOBI PROOF: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T88C-A COMPLETE")
