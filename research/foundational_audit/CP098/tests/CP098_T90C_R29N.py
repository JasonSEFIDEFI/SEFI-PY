import sympy as S
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29N ===", flush=True)

folder = Path(__file__).parent
source = (folder / "CP098_T90C_R29L.py").read_text()

# Reuse only the exact polynomial definitions from R29L.
cut = source.index('print("STAGE 2: restrict rho=2"')
prefix = source[:cut]

namespace = {"__name__": "r29n_polynomials"}
exec(compile(prefix, "R29L polynomial definitions", "exec"),
     namespace)

rho = namespace["rho"]
z = namespace["z"]
b = namespace["b"]
ND = namespace["ND"]
NT = namespace["NT"]

print("STAGE 1: reconstruct exact resultant", flush=True)

resultant = S.resultant(ND, NT, b)
assert resultant != 0

poly = S.Poly(resultant, rho, z, domain=S.ZZ)
terms = poly.terms()

assert len(terms) == 22

print("PASS: 22 exact integer terms", flush=True)

print("STAGE 2: certify even powers", flush=True)

assert all(
    erho % 2 == 0 and ez % 2 == 0
    for (erho, ez), coefficient in terms
)

print("PASS: polynomial in rho^2 and z^2", flush=True)

print("STAGE 3: certify coefficient signs", flush=True)

assert all(
    isinstance(coefficient, S.Integer)
    and coefficient > 0
    for monomial, coefficient in terms
)

constant = poly.coeff_monomial((0, 0))
print("CONSTANT COEFFICIENT =", constant, flush=True)
assert any(
    ez == 0 and coefficient > 0
    for (erho, ez), coefficient in terms
), "STOP: no positive z-independent term"

print("PASS: every coefficient strictly positive", flush=True)
print("PASS: positive z-independent term verified", flush=True)

print("STAGE 4: establish real-domain conclusion", flush=True)

print("For real rho,z with rho^2 > 31/32:")
print("  rho^2 >= 0")
print("  z^2 >= 0")
print("  resultant > 0 because rho != 0 and a positive z-independent term exists")

print("Therefore ND=NT=0 has no simultaneous")
print("real solution in the admissible domain.")

print("R29N PASS")
print("TESTED FAMILY: w=1, t=a=0")
print("D AND Ddot SIMULTANEOUS REAL ZERO: EXCLUDED")
print("GENERAL-w EXTENSION: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(perf_counter()-start, 3))
