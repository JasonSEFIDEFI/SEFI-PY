"""CP098-T56: Two-stage action-generator audit.

Test G = f(tau)*phi_T + g(tau)*chi_T.

Uses exact adapted-frame Hamiltonian bracket
identities from T47 and T49.

The coefficient functions depend only on tau.
No tertiary generator or gauge symmetry is assumed.
"""

import sympy as sp

f, g, fdot, gdot = sp.symbols(
    "f g f_dot g_dot", real=True
)

v, k = sp.symbols("v k", positive=True)

LT, LN = sp.symbols(
    "lambda_T lambda_N", real=True
)

dLT, dLN = sp.symbols(
    "delta_lambda_T delta_lambda_N",
    real=True
)

phiT, phiN, chiT, chiN, psiT = sp.symbols(
    "phi_T phi_N chi_T chi_N psi_T",
    real=True
)

print("=== CP098-T56 ===", flush=True)

# Exact adapted-frame identities:
#
# {phi_T,Hc} = -chi_T + v*k*phi_N
# {chi_T,Hc} = -psi_T + v*k*chi_N
#
# Primary brackets:
# {phi_T,phi_T} = 0
# {phi_T,phi_N} = 0
# {chi_T,phi_T} = 0
# {chi_T,phi_N} = 0

bracket_phiT_Hc = -chiT + v*k*phiN
bracket_chiT_Hc = -psiT + v*k*chiN

explicit_time = fdot*phiT + gdot*chiT

bracket_Hc = (
    f*bracket_phiT_Hc
    + g*bracket_chiT_Hc
)

# All brackets with the primary multiplier
# terms vanish for this provisional generator.
bracket_primary = sp.Integer(0)

multiplier_variation = (
    -dLT*phiT - dLN*phiN
)

bulk = sp.expand(
    explicit_time
    + bracket_Hc
    + bracket_primary
    + multiplier_variation
)

expected = (
    (fdot-dLT)*phiT
    + (f*v*k-dLN)*phiN
    + (gdot-f)*chiT
    + g*v*k*chiN
    - g*psiT
)

assert sp.expand(bulk-expected) == 0

print("Full two-stage bulk variation:")
sp.pprint(bulk)

print("PASS: two-stage action identity",
      flush=True)

# The secondary chi_T coefficient must vanish.
chiT_coefficient = sp.diff(bulk,chiT)

assert sp.factor(
    chiT_coefficient-(gdot-f)
) == 0

print("chi_T coefficient =", chiT_coefficient,
      flush=True)

# Necessary chain condition:
# f = gdot.

chain_condition = {f:gdot}

after_chain = sp.factor(
    bulk.subs(chain_condition)
)

assert sp.diff(after_chain,chiT) == 0

print("PASS: secondary obstruction canceled",
      flush=True)

# Primary multiplier variations can cancel
# the remaining primary terms.
primary_cancellation = {
    dLT:fdot,
    dLN:gdot*v*k
}

remaining = sp.factor(
    after_chain.subs(primary_cancellation)
)

expected_remaining = (
    g*v*k*chiN-g*psiT
)

assert sp.expand(
    remaining-expected_remaining
) == 0

print("After chain and primary cancellation:")
sp.pprint(remaining)

print("PASS: remaining secondary-normal "
      "and tertiary terms identified",
      flush=True)

# On the primary + secondary surface,
# chi_N=0, but psi_T is not yet imposed.
primary_secondary = {
    phiT:0,
    phiN:0,
    chiT:0,
    chiN:0
}

reduced = sp.factor(
    remaining.subs(primary_secondary)
)

assert sp.factor(reduced+g*psiT) == 0

print("Primary-secondary residual =", reduced,
      flush=True)

# For arbitrary nonzero g, this residual
# cannot be canceled by varying only
# the primary multipliers.
assert sp.diff(reduced,psiT) == -g

print("PASS: tertiary obstruction persists",
      flush=True)

# On the complete five-constraint surface
# the variation vanishes weakly.
complete = dict(primary_secondary)
complete[psiT] = 0

assert sp.factor(
    remaining.subs(complete)
) == 0

print("PASS: weak vanishing distinguished",
      flush=True)

# Important: the normal secondary term also
# remains off shell. A complete generator
# may need normal-sector corrections.
assert sp.diff(remaining,chiN) == g*v*k

print("PASS: normal-sector correction needed",
      flush=True)

print("TWO-STAGE GENERATOR: INSUFFICIENT OFF SHELL")
print("NECESSARY RELATION: f = g_dot")
print("TERTIARY TERM: -g*psi_T")
print("NORMAL SECONDARY TERM: g*v*k*chi_N")
print("FULL GENERATOR: OPEN")
print("ACTION INVARIANCE: NOT ESTABLISHED")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
