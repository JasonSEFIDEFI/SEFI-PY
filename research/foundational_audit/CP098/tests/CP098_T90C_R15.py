import sympy as sp

print("=== CP098-T90C-R15 ===")

eta = sp.diag(-1, 1, 1, 1)

def dot(a, b):
    return (a.T * eta * b)[0]

def simplify_matrix(M):
    return M.applyfunc(sp.simplify)

A, B, C = sp.symbols("A B C", real=True, nonzero=True)

U = sp.symbols("u0:4", real=True)
Q = sp.symbols("q0:4", real=True)
J = sp.symbols("j0:4", real=True)

u = sp.Matrix(U)
q = sp.Matrix(Q)
j = sp.Matrix(J)

v = sp.sqrt(-dot(u, u))
T = u/v

w = q + dot(T, q)*T
L = sp.sqrt(dot(w, w))
k = L/v**2
N = w/L

projector = (
    sp.eye(4)
    + T*(T.T*eta)
    - N*(N.T*eta)
)

jt = projector*j

# Candidate action terms derived from geometry.
Lcurv = -v*A*k**2
Ltors = -C*dot(jt,jt)/(v**5*k**2)

# B includes the mass and internal Hamiltonian.
Lcandidate = -v*(B+A*k**2) + Ltors

# Curvature derivative.
pcurv = sp.Matrix([
    sp.diff(Lcurv, qi) for qi in Q
])

expected_pcurv = -2*A*k/v * eta*N

# Highest-derivative momentum.
P = sp.Matrix([
    sp.diff(Ltors, ji) for ji in J
])

P2 = sp.simplify(dot(eta*P, eta*P))

# Highest-derivative Legendre expression.
Hhighest = (P.T*j)[0] - Lcandidate

expected_H = (
    v*(B+A*k**2)
    -v**5*k**2*P2/(4*C)
)

jets = [
    (
        [5,3,0,0],
        [1,2,2,1],
        [4,1,3,6]
    ),
    (
        [13,5,12,0],
        [2,3,1,4],
        [1,2,5,7]
    )
]

for index, (uv,qv,jv) in enumerate(jets, 1):
    subs = dict(zip(U, map(sp.Integer, uv)))
    subs.update(dict(zip(Q, map(sp.Integer, qv))))
    subs.update(dict(zip(J, map(sp.Integer, jv))))

    v0 = sp.simplify(v.subs(subs))
    k0 = sp.simplify(k.subs(subs))

    assert v0.is_positive
    assert k0.is_positive

    curvature_residual = simplify_matrix(
        (pcurv-expected_pcurv).subs(subs)
    )

    assert curvature_residual == sp.zeros(4,1)

    legendre_residual = sp.simplify(
        (Hhighest-expected_H).subs(subs)
    )

    assert legendre_residual == 0

    print()
    print("JET", index)
    print("v =", v0)
    print("kappa =", k0)
    print("CURVATURE RESIDUAL =", curvature_residual)
    print("LEGENDRE RESIDUAL =", legendre_residual)
    print("PASS")

print()
print("R15 COMPLETE")
print("CURVATURE DERIVATIVE: VERIFIED AT TWO JETS")
print("TORSION LEGENDRE IDENTITY: VERIFIED AT TWO JETS")
print("FULL ACTION EQUIVALENCE: NOT ESTABLISHED")
