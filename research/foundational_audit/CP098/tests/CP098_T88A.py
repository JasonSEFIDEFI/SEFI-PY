"""CP098-T88A: mixed-second-derivative prerequisite.

Exact non-adapted points.
Checks independent orders of differentiation.
Does not establish the full constraint Jacobi identity.
"""
import sympy as sp

print("=== CP098-T88A ===", flush=True)

eta = sp.diag(-1, 1, 1, 1)
e, f = sp.symbols("e f", real=True)

def dot(x, y):
    return (x.T * eta * y)[0]

def geometry(u, q):
    v = sp.sqrt(-dot(u, u))
    T = u/v
    w = q + T*dot(T, q)
    L = sp.sqrt(dot(w, w))
    k = L/v**2
    return {"v": v, "L": L, "kappa": k}

cases = [
    (sp.Matrix([5,3,0,0]), sp.Matrix([3,5,0,0])),
    (sp.Matrix([10,2,3,1]), sp.Matrix([4,1,5,2])),
]

directions = [
    (kind, i)
    for kind in ("u", "q")
    for i in range(4)
]

def perturb(u, q, direction, amount):
    kind, i = direction
    h = sp.eye(4)[:,i]*amount
    if kind == "u":
        return u+h, q
    return u, q+h

total = 0

for case_id, (u0,q0) in enumerate(cases,1):
    for ia, da in enumerate(directions):
        for db in directions[ia+1:]:
            u1,q1 = perturb(u0,q0,da,e)
            u2,q2 = perturb(u1,q1,db,f)

            values = geometry(u2,q2)

            for name, expression in values.items():
                first = sp.diff(
                    sp.diff(expression,e),f
                ).subs({e:0,f:0})

                second = sp.diff(
                    sp.diff(expression,f),e
                ).subs({e:0,f:0})

                residual = sp.simplify(first-second)

                assert residual == 0, (
                    f"STOP: case={case_id}, "
                    f"directions={da},{db}, "
                    f"quantity={name}, "
                    f"residual={residual}"
                )
                total += 1

    print(
        f"PASS: case {case_id}, mixed derivatives",
        flush=True
    )

print("MIXED_DERIVATIVE_IDENTITIES =", total)
assert total == 168

print("PASS: geometric mixed-derivative symmetry")
print("FULL CONSTRAINT JACOBI AUDIT: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T88A COMPLETE")
