import math
import numpy as np

print("Checkpoint 022: Singular-value geometry")

s = 0.01

def cycle(u,v):
    return math.exp(-(u+v))*np.array([
        [1.0,10*u],
        [10*v,1.0+100*u*v]
    ])

def spectral_radius(M):
    return max(abs(np.linalg.eigvals(M)))

def largest_singular(M):
    return np.linalg.svd(M,compute_uv=False)[0]

results = []

for L in (1.0,2.0,3.0,4.0):
    B = cycle(s,L)
    C = cycle(L,s)
    T = C @ B

    assert np.allclose(C,B.T,rtol=1e-13,atol=1e-13)
    assert np.allclose(T,B.T@B,rtol=1e-13,atol=1e-13)

    rho_B = spectral_radius(B)
    sigma_B = largest_singular(B)
    rho_T = spectral_radius(T)

    assert np.isclose(
        rho_T,sigma_B**2,rtol=1e-12
    )

    assert rho_B < 1

    results.append((L,rho_B,sigma_B,rho_T))

    print(
        f"L={L:.1f} "
        f"rho_B={rho_B:.9f} "
        f"sigma_B={sigma_B:.9f} "
        f"rho_CB={rho_T:.9f}"
    )

print("PASS: exact transpose and singular relations")

# The first three alternating pairs are unstable.
for L,rho,sigma,rho_T in results[:3]:
    assert sigma > 1
    assert rho_T > 1

# The fourth alternating pair is contractive
# in the Euclidean norm.
L,rho,sigma,rho_T = results[3]

assert sigma < 1
assert rho_T < 1

print("PASS: stable/unstable transition bracketed")

# At L=4, B and C share the Euclidean
# quadratic Lyapunov function V(x)=x.T@x.
B = cycle(s,4.0)
C = B.T
I = np.eye(2)

for M in (B,C):
    difference = M.T @ M-I
    assert max(np.linalg.eigvalsh(difference)) < 0

print("PASS: common Euclidean Lyapunov certificate")

# At L=1, no common strict quadratic
# Lyapunov certificate can exist, because
# the alternating product has rho > 1.
B = cycle(s,1.0)
C = B.T
assert spectral_radius(C@B) > 1

print("PASS: common-certificate obstruction")
print("PASS: Checkpoint 022 verified")
