import sympy as sp

eta = sp.diag(-1, 1, 1, 1)

def dot(a, b):
    return sp.simplify((a.T * eta * b)[0])

def check(u, q, j):
    v = sp.sqrt(-dot(u, u))
    T = u / v

    w = q + dot(T, q) * T
    L = sp.sqrt(dot(w, w))
    N = w / L

    kappa = sp.simplify(L / v**2)

    def transverse(X):
        return sp.simplify(
            X + dot(T, X)*T - dot(N, X)*N
        )

    # Differentiate N(u,q) along the jet:
    # du/dlambda=q, dq/dlambda=j.
    U = sp.symbols("u0:4")
    Q = sp.symbols("q0:4")
    um = sp.Matrix(U)
    qm = sp.Matrix(Q)

    vm = sp.sqrt(-dot(um, um))
    Tm = um/vm
    wm = qm + dot(Tm, qm)*Tm
    Lm = sp.sqrt(dot(wm, wm))
    Nm = wm/Lm

    dN = sp.zeros(4, 1)
    for a in range(4):
        dN += Nm.diff(U[a])*q[a]
        dN += Nm.diff(Q[a])*j[a]

    subs = dict(zip(U, list(u)))
    subs.update(dict(zip(Q, list(q))))

    dN = sp.simplify(dN.subs(subs))

    lhs = sp.simplify(
        dot(transverse(dN/v), transverse(dN/v))
    )

    rhs = sp.simplify(
        dot(transverse(j), transverse(j))
        / (v**6 * kappa**2)
    )

    residual = sp.simplify(lhs-rhs)

    print("v =", v)
    print("kappa =", kappa)
    print("torsion squared =", lhs)
    print("formula =", rhs)
    print("residual =", residual)

    assert residual == 0
    assert sp.simplify(dot(T,T)+1) == 0
    assert sp.simplify(dot(N,N)-1) == 0
    assert sp.simplify(dot(T,N)) == 0

    return lhs

print("=== CP098-T90C-R13 ===")

jets = [
    (
        sp.Matrix([3,0,0,0]),
        sp.Matrix([1,18,0,0]),
        sp.Matrix([2,3,5,7])
    ),
    (
        sp.Matrix([5,3,0,0]),
        sp.Matrix([1,2,2,1]),
        sp.Matrix([4,1,3,6])
    ),
]

for index, (u,q,j) in enumerate(jets, 1):
    print()
    print("JET", index)
    result = check(u,q,j)
    print("PASS: exact torsion identity")

print()
print("R13 COMPLETE")
print("LIMIT: kinematic identity only")
print("DEFI identification: unresolved")
print("LEGENDRE DERIVATION: NOT YET PERFORMED")
