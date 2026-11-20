"""CP098-T88B-R3: full canonical Jacobi audit.

Triple: (phi_T, phi_N, chi_N).
Two exact non-adapted points, unrestricted momenta.
All canonical (u,p1) and (q,P) contributions retained.
"""
import sympy as sp

print("=== CP098-T88B-R3 ===", flush=True)

eta = sp.diag(-1, 1, 1, 1)
e = sp.symbols("e", real=True)
basis = sp.eye(4)
A = sp.Rational(3, 2)

def dot(x, y):
    return (x.T*eta*y)[0]

def pair(p, x):
    return (p.T*x)[0]

def geom(u, q):
    v = sp.sqrt(-dot(u, u))
    T = u/v
    w = q+T*dot(T, q)
    L = sp.sqrt(dot(w, w))
    N = w/L
    return v,T,N,L

def phiT(u,q,P,p1):
    return pair(P,geom(u,q)[1])

def phiN(u,q,P,p1):
    return pair(P,geom(u,q)[2])

def chiN(u,q,P,p1):
    v,T,N,L = geom(u,q)
    return pair(p1,N)+2*A*L/v**3

def diffdir(fn,args,which,h):
    vals = list(args)
    vals[which] = vals[which]+e*h
    return sp.diff(fn(*vals),e).subs(e,0)

def pb(f,g,args):
    total = sp.S.Zero
    for i in range(4):
        h = basis[:,i]
        fu = diffdir(f,args,0,h)
        fq = diffdir(f,args,1,h)
        fP = diffdir(f,args,2,h)
        fj = diffdir(f,args,3,h)
        gu = diffdir(g,args,0,h)
        gq = diffdir(g,args,1,h)
        gP = diffdir(g,args,2,h)
        gj = diffdir(g,args,3,h)

        total += fu*gj-fj*gu+fq*gP-fP*gq

    return sp.simplify(total)

def b_Nchi(u,q,P,p1):
    v,T,N,L = geom(u,q)
    return -2*A/v**3+phiT(u,q,P,p1)/v

def b_chiT(u,q,P,p1):
    v,T,N,L = geom(u,q)
    return -phiN(u,q,P,p1)/v

def zero(u,q,P,p1):
    return sp.S.Zero

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

P = sp.Matrix([2,3,5,7])
p1 = sp.Matrix([11,13,17,19])

for case_id,(u,q) in enumerate(cases,1):
    args = (u,q,P,p1)

    inner = [
        ("PHIT_PHIN",phiT,phiN,zero),
        ("PHIN_CHIN",phiN,chiN,b_Nchi),
        ("CHIN_PHIT",chiN,phiT,b_chiT)
    ]

    for label,f,g,expected in inner:
        residual = sp.simplify(
            pb(f,g,args)-expected(*args)
        )
        print(
            f"CASE_{case_id}_{label}_RESIDUAL =",
            residual,flush=True
        )
        assert residual == 0, (
            f"STOP: {label} failed"
        )

    # Actual outer canonical brackets, not a
    # preassigned zero or a Jacobi assumption.
    term1 = pb(phiT,b_Nchi,args)
    term2 = pb(phiN,b_chiT,args)
    term3 = pb(chiN,zero,args)

    for i,term in enumerate(
        (term1,term2,term3),1
    ):
        print(
            f"CASE_{case_id}_OUTER_{i} =",
            term,flush=True
        )

    jacobi = sp.simplify(
        term1+term2+term3
    )

    print(
        f"CASE_{case_id}_JACOBI_RESIDUAL =",
        jacobi,flush=True
    )
    assert jacobi == 0, (
        f"STOP: Jacobi mismatch in case {case_id}"
    )

    print(
        f"PASS: full canonical Jacobi case {case_id}",
        flush=True
    )

print("PASS: corrected full inner-bracket checks")
print("PASS: explicitly evaluated outer brackets")
print("ARBITRARY-POINT JACOBI PROOF: OPEN")
print("PSI_T JACOBI TRIPLES: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T88B-R3 COMPLETE")
