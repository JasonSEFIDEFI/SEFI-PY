#!/usr/bin/env python3
"""
GWFM / SEFI — CHECKPOINT 091

DYNAMICAL PRESERVATION OF THE CLASSICAL
MULTIPLICATION OBSERVABLE ALGEBRA

Classical free relativistic flow:
    F_t(x,p) = (x + v(p)*t, p)

Koopman action:
    U_t psi = psi composed with F_{-t}

Heisenberg action:
    U_{-t} M_f U_t = M_{f composed with F_t}

The analytic identity applies to bounded
measurable multiplication observables.

Symbolic tests use smooth representatives.

No quantum observable or measurement
postulate is introduced.
"""

import sympy as sp

I = sp.I
x, p, t, s = sp.symbols(
    "x p t s", real=True
)
m = sp.symbols("m", positive=True)
a, b, c = sp.symbols("a b c", real=True)

v = p/sp.sqrt(p**2+m**2)

psi = sp.exp(-x**2-p**2)
f = x**2+2*x*p+p**2
g = 3*x-p**2
chi = a*x*p+b*x+c*p**2

print("="*64)
print("GWFM / SEFI — CHECKPOINT 091")
print("DYNAMICAL PRESERVATION OF CLASSICAL ALGEBRA")
print("="*64)

def eq(left,right):
    return sp.simplify(
        sp.expand(left-right)
    ) == 0

def U(expr,time):
    return expr.subs(
        x,x-v*time,simultaneous=True
    )

def evolved_observable(expr,time):
    return expr.subs(
        x,x+v*time,simultaneous=True
    )

def M(obs,state):
    return obs*state

def phase(state,phase_function):
    return sp.exp(I*phase_function)*state

def density(state):
    return sp.simplify(
        state*sp.conjugate(state)
    )

# 1. The flow has a group law.
assert eq(
    U(U(x,s),t),
    U(x,s+t)
)
assert eq(
    U(U(p,s),t),
    p
)

print("PASS: free relativistic flow")
print("      satisfies the group law")

# 2. Inverse flow.
assert eq(
    U(U(x,t),-t),x
)
assert eq(
    U(U(p,t),-t),p
)

print("PASS: inverse flow exists")

# 3. Phase-space Jacobian.
F = sp.Matrix([x+v*t,p])
J = F.jacobian([x,p])

assert eq(J.det(),1)

print("PASS: flow preserves")
print("      phase-space volume")

# 4. Exact Heisenberg evolution
# of a multiplication observable.
#
# Use an independent test function
# so the result is not dependent
# on Gaussian simplification.

test = x**3+x*p+2*p**2

left = U(
    M(f,U(test,t)),
    -t
)

right = M(
    evolved_observable(f,t),
    test
)

assert eq(left,right)

print("PASS: conjugation of a")
print("      multiplication operator")
print("      remains multiplication")

# 5. The identity holds for several
# distinct polynomial observables.
for obs in (x,p,x*p,x**2+p**2,f,g):
    left = U(
        M(obs,U(test,t)),
        -t
    )
    right = M(
        evolved_observable(obs,t),
        test
    )
    assert eq(left,right)

print("PASS: covariance verified")
print("      for multiple observables")

# 6. The observable evolution
# preserves multiplication.
assert eq(
    evolved_observable(f*g,t),
    evolved_observable(f,t)
    *evolved_observable(g,t)
)

print("PASS: evolution preserves")
print("      observable products")

# 7. It also preserves adjoints.
complex_obs = f+I*g

assert eq(
    evolved_observable(
        sp.conjugate(complex_obs),t
    ),
    sp.conjugate(
        evolved_observable(complex_obs,t)
    )
)

print("PASS: evolution preserves")
print("      observable adjoints")

# 8. Classical observable evolution
# is commutative.
ft = evolved_observable(f,t)
gt = evolved_observable(g,t)

assert eq(ft*gt,gt*ft)

print("PASS: evolved classical algebra")
print("      remains commutative")

# 9. Local phase transformations
# remain multiplication operators.
phase_t = evolved_observable(chi,t)

assert eq(
    evolved_observable(
        sp.exp(I*chi),t
    ),
    sp.exp(I*phase_t)
)

print("PASS: local phase operators")
print("      remain local phase operators")

# 10. Exact covariance of a phase
# under the Koopman representation.
#
# U_-t W_chi U_t = W_(chi composed F_t).

left = U(
    phase(U(test,t),chi),
    -t
)

right = phase(test,phase_t)

assert eq(left,right)

print("PASS: phase transformations")
print("      transform covariantly")

# 11. Density invariance under
# real local phase multiplication.
assert eq(
    density(phase(psi,chi)),
    density(psi)
)

print("PASS: local phases preserve")
print("      the classical density")

# 12. Density equivalence is
# preserved by the classical flow.
psi_t = U(psi,t)
phi_t = U(phase(psi,chi),t)

assert eq(
    density(psi_t),
    density(phi_t)
)

print("PASS: density equivalence")
print("      preserved by evolution")

# 13. Liouville transport.
rho_t = density(psi_t)

assert eq(
    sp.diff(rho_t,t)
    +v*sp.diff(rho_t,x),
    0
)

print("PASS: density obeys")
print("      Liouville transport")

# 14. The infinitesimal generator
# acts as the classical derivation.
assert eq(
    sp.diff(ft,t).subs(t,0),
    v*sp.diff(f,x)
)

print("PASS: infinitesimal observable")
print("      evolution is classical")

# 15. Coordinate dynamics.
assert eq(
    evolved_observable(x,t),
    x+v*t
)

assert eq(
    evolved_observable(p,t),
    p
)

print("PASS: classical position")
print("      and momentum laws recovered")

# 16. No special phase scale
# enters the classical covariance.
lam = sp.symbols(
    "lam",real=True,nonzero=True
)

scaled_phase = chi/lam

assert eq(
    evolved_observable(
        scaled_phase,t
    ),
    phase_t/lam
)

print("PASS: classical covariance")
print("      selects no phase scale")

print("\n"+"="*64)
print("CHECKPOINT 091: ALL TESTS PASSED")
print("="*64)

print("ESTABLISHED MATHEMATICALLY:")
print("- The free classical flow is invertible")
print("  and preserves phase-space volume.")
print("- Koopman conjugation preserves the")
print("  multiplication observable algebra.")
print("- The evolved algebra remains")
print("  commutative and closed under adjoints.")
print("- Local phase multiplication transforms")
print("  covariantly under classical evolution.")
print("- Density equivalence is preserved.")
print("- The induced density evolution")
print("  satisfies the Liouville equation.")
print("- No preferred phase scale is selected.")

print("\nANALYTIC SCOPE:")
print("- For bounded measurable f, the exact")
print("  Koopman conjugation identity follows")
print("  from composition with the invertible")
print("  measure-preserving classical flow.")
print("- Polynomial examples are generally")
print("  unbounded; their symbolic identities")
print("  hold on suitable common test domains.")
print("- The tests do not prove a new physical")
print("  observable-selection principle.")

print("\nNOT ESTABLISHED:")
print("- A physically selected noncommutative")
print("  observable algebra.")
print("- Quantum measurement postulates.")
print("- A preferred nonzero action scale.")
print("- Quantum dynamics emerging from GWFM.")
print("- A derivation of quantum mechanics.")
