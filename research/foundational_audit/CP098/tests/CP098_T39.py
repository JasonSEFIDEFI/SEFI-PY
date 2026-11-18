"""CP098-T39: Regular-sector closure assembly.

Assembles previously tested bracket identities.
This is a consistency certificate, not a new
independent canonical differentiation.

Regular sector: v>0, kap>0, A!=0, C!=0.
"""

import sympy as sp

v, kap = sp.symbols("v kap", positive=True)
A, C = sp.symbols("A C", nonzero=True, real=True)
DA = sp.symbols("DA", real=True)
FN_ext = sp.symbols("FN_ext", real=True)
LT = sp.symbols("lambda_T", real=True)

FN = FN_ext + 2*kap*DA
LN = -v**3*FN/(2*A)

# Each expression below is supported by an
# earlier independent canonical calculation.
#
# T38R: primary preservation and primary brackets.
dot_phi_T = sp.Integer(0)
dot_phi_N = sp.Integer(0)

# T37: tangential secondary preservation.
dot_chi_T = sp.Integer(0)

# T34: normal secondary preservation.
dot_chi_N = FN + LN*(2*A/v**3)

# T35: external tangential tertiary bracket.
psi_H_external = -v*kap*FN_ext

# T36: internal tangential tertiary bracket.
psi_H_internal = -2*v*kap**2*DA

# T35: multiplier brackets.
psi_phi_T = sp.Integer(0)
psi_phi_N = -2*A*kap/v**2

dot_psi_T = (
    psi_H_external
    + psi_H_internal
    + LT*psi_phi_T
    + LN*psi_phi_N
)

equations = {
    "dot(phi_T)": dot_phi_T,
    "dot(phi_N)": dot_phi_N,
    "dot(chi_T)": dot_chi_T,
    "dot(chi_N)": dot_chi_N,
    "dot(psi_T)": dot_psi_T,
}

print("=== CP098-T39 ===", flush=True)
print("Regular sector: v>0, kap>0, A!=0, C!=0")
print("Normal multiplier =", sp.factor(LN))

for name, expr in equations.items():
    residual = sp.factor(expr)
    print(f"{name} = {residual}", flush=True)
    assert residual == 0, name

print("PASS: all five assembled preservation equations")
print("PASS: one consistent normal multiplier")
print("lambda_T: undetermined by these equations")
print("REGULAR-SECTOR CLOSURE: SUPPORTED")
print("INDEPENDENT DERIVATION: T34-T38R")
print("SINGULAR-SECTOR CLOSURE: OPEN")
print("CONSTRAINT CLASSIFICATION: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
