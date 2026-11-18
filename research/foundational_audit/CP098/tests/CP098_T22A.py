"""CP098-T22A: Independent transverse-factor audit.

Verifies the v^2 coefficient in the q,P contribution.
Does not establish general Dirac closure.
"""

import sympy as sp

v, k, C = sp.symbols("v k C", positive=True)
r, s, x, y = sp.symbols("r s x y", real=True)

# At the adapted regular point:
# w^2 = v^4*k^2.
w2 = v**4*k**2

dpsi_dP2 = -w2*r/(2*C)
dpsi_dP3 = -w2*s/(2*C)

dchi_dq2 = x/(v**2*k)
dchi_dq3 = y/(v**2*k)

term2 = sp.factor(-dpsi_dP2*dchi_dq2)
term3 = sp.factor(-dpsi_dP3*dchi_dq3)

actual = sp.factor(term2 + term3)
expected = v**2*k*(r*x+s*y)/(2*C)

print("=== CP098-T22A ===")
print("Component 2 =", term2)
print("Component 3 =", term3)
print("Combined =", actual)
print("Expected corrected =", expected)

assert sp.simplify(actual-expected) == 0

# Check that the original candidate really differs
# away from special values.
old = v*k*(r*x+s*y)/(2*C)
difference = sp.factor(actual-old)

print("Difference from old candidate =", difference)

assert difference != 0

print("PASS: corrected v^2 coefficient")
print("OLD DIAGNOSTIC CANDIDATE: REJECTED")
print("GENERAL DIRAC CLOSURE: UNRESOLVED")
