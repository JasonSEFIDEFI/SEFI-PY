"""CP098-T86: Independent non-adapted canonical bracket checks.

Two exact geometric test points.
All canonical momentum components unrestricted.
No momentum constraints imposed.
Not a general-point symbolic proof.
"""
import sympy as sp

print("=== CP098-T86 ===", flush=True)

eta = sp.diag(-1,1,1,1)
e = sp.symbols("e", real=True)

A = sp.Rational(3,2)
C = sp.Rational(5,3)
B = sp.Rational(7,4)
DA = sp.Rational(2,5)

P = sp.Matrix([2,3,5,7])
p1 = sp.Matrix([11,13,17,19])
p0 = sp.Matrix([23,29,31,37])

basis = sp.eye(4)

def dot(X,Y):
    return (X.T*eta*Y)[0]

def pair(p,X):
    return (p.T*X)[0]

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

def chi(u,q):
    v,T,N,L,k = geom(u,q)
    return pair(p1,N)+2*A*k/v

def psi(u,q):
    v,T,N,L,k = geom(u,q)
    return pair(p0,T)+B-A*k**2-v**4*k**2*dot(P,P)/(4*C)

def compact(u,q):
    v,T,N,L,k = geom(u,q)
    tq = dot(T,q)
    return (
        -L/(2*C)*(
            covdot(p1,P)
            +pair(p1,T)*covdot(T,P)
            -pair(p1,N)*covdot(N,P)
        )
        +L*pair(p1,T)*covdot(P,T)/C
        -A*L**2*covdot(N,P)/(C*v**3)
        -pair(p0,N)/v
        +2*A*L*tq/v**5
        +L*tq*dot(P,P)/(2*C*v)
        +2*L*DA/v**3
    )

def directional(f,u,q,which,i):
    if which=="u":
        expression=f(u+e*basis[:,i],q)
    else:
        expression=f(u,q+e*basis[:,i])
    return sp.diff(expression,e).subs(e,0)

cases = [
    (
        sp.Matrix([5,3,0,0]),
        sp.Matrix([3,5,0,0])
    ),
    (
        sp.Matrix([10,2,3,1]),
        sp.Matrix([4,1,5,2])
    )
]

for index,(u,q) in enumerate(cases,1):
    v,T,N,L,k = geom(u,q)

    assert sp.simplify(-dot(u,u))>0
    assert sp.simplify(dot(q+T*dot(T,q),
                           q+T*dot(T,q)))>0

    E_external = 0

    for i in range(4):
        chi_q = directional(chi,u,q,"q",i)
        psi_P = -v**4*k**2*(eta*P)[i]/(2*C)
        psi_u = directional(psi,u,q,"u",i)

        E_external += chi_q*psi_P-N[i]*psi_u

    E_direct = E_external+2*k*DA/v
    E_compact = compact(u,q)

    residual = sp.simplify(
        sp.cancel(sp.together(E_direct-E_compact))
    )

    print(f"CASE_{index}_RESIDUAL =",residual,
          flush=True)

    assert residual==0, (
        f"STOP: canonical bracket mismatch in case {index}"
    )

    print(f"PASS: non-adapted case {index}",
          flush=True)

print("EXACT GENERAL-POINT SPOT CHECKS: VERIFIED")
print("ARBITRARY-POINT SYMBOLIC PROOF: OPEN")
print("JACOBI CROSS-CHECK: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T86 COMPLETE")
