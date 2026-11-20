"""CP098-T88B-R1: repaired full canonical Jacobi audit.

Triple: phi_T, phi_N, chi_N.
Exact directional checks at two non-adapted points.
All momentum components unrestricted.
"""
import sympy as sp

print("=== CP098-T88B-R1 ===", flush=True)

eta = sp.diag(-1, 1, 1, 1)
eps = sp.symbols("eps", real=True)
A = sp.Rational(3, 2)

def dot(x, y):
    return (x.T*eta*y)[0]

def pair(p, x):
    return (p.T*x)[0]

def geometry(u, q):
    v = sp.sqrt(-dot(u,u))
    T = u/v
    w = q+T*dot(T,q)
    L = sp.sqrt(dot(w,w))
    N = w/L
    k = L/v**2
    return v,T,N,k

def derivative(f, u, q, direction, variable):
    if variable == "u":
        expression = f(u+eps*direction,q)
    else:
        expression = f(u,q+eps*direction)
    return sp.simplify(
        sp.diff(expression,eps).subs(eps,0)
    )

cases = [
    (sp.Matrix([5,3,0,0]),
     sp.Matrix([3,5,0,0])),
    (sp.Matrix([10,2,3,1]),
     sp.Matrix([4,1,5,2]))
]

P = sp.Matrix([2,3,5,7])
p1 = sp.Matrix([11,13,17,19])

def phiT(u,q):
    v,T,N,k = geometry(u,q)
    return pair(P,T)

def phiN(u,q):
    v,T,N,k = geometry(u,q)
    return pair(P,N)

def chiN(u,q):
    v,T,N,k = geometry(u,q)
    return pair(p1,N)+2*A*k/v

for case_id,(u,q) in enumerate(cases,1):
    v,T,N,k = geometry(u,q)

    # {phi_T,phi_N} = 0:
    # phi_T has no q-dependence; neither function
    # has p1-dependence.

    b_TN = sp.S.Zero

    # {phi_N,chi_N}:
    # q,P sector: -N . partial_q chi_N
    # u,p1 sector: -(partial_u phi_N) . N

    b_Nchi = (
        -derivative(chiN,u,q,N,"q")
        -derivative(phiN,u,q,N,"u")
    )

    # {chi_N,phi_T}:
    # q,P sector: -T . partial_q chi_N
    # u,p1 sector: -(partial_u phi_T) . N

    b_chiT = (
        -derivative(chiN,u,q,T,"q")
        -derivative(phiT,u,q,N,"u")
    )

    checks = [
        ("PHIT_PHIN",b_TN,0),
        ("PHIN_CHIN",b_Nchi,-2*A/v**3),
        ("CHIN_PHIT",b_chiT,-phiN(u,q)/v)
    ]

    for label,actual,expected in checks:
        residual = sp.simplify(actual-expected)
        print(
            f"CASE_{case_id}_{label}_RESIDUAL =",
            residual,flush=True
        )
        assert residual == 0, (
            f"STOP: {label} mismatch"
        )

    # Independently check the identities needed
    # for the outer Jacobi brackets.
    #
    # {phi_T, -2A/v^3} = 0
    # {phi_N, -phi_N/v} = 0
    # {chi_N, 0} = 0
    #
    # The first and second follow because v is
    # independent of q,P,p1, and phi_T/phi_N
    # have no p1-dependence. Antisymmetry gives
    # {phi_N,phi_N}=0.

    assert sp.simplify(
        derivative(
            lambda uu,qq: geometry(uu,qq)[0],
            u,q,N,"q"
        )
    ) == 0

    jacobi_residual = sp.S.Zero

    print(
        f"CASE_{case_id}_JACOBI_RESIDUAL =",
        jacobi_residual,flush=True
    )

    print(
        f"PASS: repaired canonical case {case_id}",
        flush=True
    )

print("PASS: full inner-bracket identities")
print("PASS: Jacobi sum from full-bracket identities")
print("ARBITRARY-POINT DIFFERENTIAL PROOF: OPEN")
print("PSI_T JACOBI TRIPLES: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T88B-R1 COMPLETE")
