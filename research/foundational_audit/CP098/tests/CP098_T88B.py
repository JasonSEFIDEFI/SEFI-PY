"""CP098-T88B: first constraint-specific Jacobi audit.

Triple: (phi_T, phi_N, chi_N).
Exact non-adapted regular points.
Unrestricted momenta.
"""
import sympy as sp

print("=== CP098-T88B ===", flush=True)

eta = sp.diag(-1, 1, 1, 1)

def dot(x, y):
    return (x.T*eta*y)[0]

def pair(p, x):
    return (p.T*x)[0]

def geom(u, q):
    v = sp.sqrt(-dot(u, u))
    T = u/v
    w = q+T*dot(T,q)
    L = sp.sqrt(dot(w,w))
    N = w/L
    k = L/v**2
    return v,T,N,k

def canonical_qP_bracket(f, g, q, P):
    return sum(
        sp.diff(f,q[i])*sp.diff(g,P[i])
        -sp.diff(f,P[i])*sp.diff(g,q[i])
        for i in range(4)
    )

u_values = [
    [5,3,0,0],
    [10,2,3,1]
]
q_values = [
    [3,5,0,0],
    [4,1,5,2]
]

P_values = [2,3,5,7]
p1_values = [11,13,17,19]

A = sp.Rational(3,2)

# Differentiate before substituting the geometric point.
q = sp.Matrix(sp.symbols("q0:4", real=True))
P = sp.Matrix(sp.symbols("P0:4", real=True))
p1 = sp.Matrix(sp.symbols("j0:4", real=True))

for case_id,(uv,qv) in enumerate(
    zip(u_values,q_values),1
):
    u = sp.Matrix(uv)
    v,T,N,k = geom(u,q)

    phiT = pair(P,T)
    phiN = pair(P,N)
    chiN = pair(p1,N)+2*A*k/v

    b_TN = canonical_qP_bracket(
        phiT,phiN,q,P
    )
    b_Nchi = canonical_qP_bracket(
        phiN,chiN,q,P
    )
    b_chiT = canonical_qP_bracket(
        chiN,phiT,q,P
    )

    subs = dict(zip(q,qv))
    subs.update(dict(zip(P,P_values)))
    subs.update(dict(zip(p1,p1_values)))

    expected_Nchi = -2*A/v**3

    checks = [
        ("PHIT_PHIN",b_TN,sp.S.Zero),
        ("PHIN_CHIN",b_Nchi,expected_Nchi),
        ("CHIN_PHIT",b_chiT,sp.S.Zero),
    ]

    for label,actual,expected in checks:
        residual = sp.simplify(
            (actual-expected).subs(subs)
        )
        print(
            f"CASE_{case_id}_{label}_RESIDUAL =",
            residual,flush=True
        )
        assert residual == 0, (
            f"STOP: {label} mismatch"
        )

    jacobi = (
        canonical_qP_bracket(
            phiT,b_Nchi,q,P
        )
        +canonical_qP_bracket(
            phiN,b_chiT,q,P
        )
        +canonical_qP_bracket(
            chiN,b_TN,q,P
        )
    )

    residual = sp.simplify(jacobi.subs(subs))

    print(
        f"CASE_{case_id}_JACOBI_RESIDUAL =",
        residual,flush=True
    )

    assert residual == 0, (
        f"STOP: Jacobi identity failed "
        f"in case {case_id}"
    )

    print(
        f"PASS: constraint Jacobi case {case_id}",
        flush=True
    )

print("PASS: first constraint-specific Jacobi audit")
print("GENERAL-POINT JACOBI PROOF: OPEN")
print("JACOBI TRIPLES INVOLVING PSI_T: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T88B COMPLETE")
