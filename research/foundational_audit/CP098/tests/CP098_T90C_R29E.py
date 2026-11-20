import sympy as S

print("=== CP098-T90C-R29E ===", flush=True)

w, rho, z, b = S.symbols(
    "w rho z b", real=True
)

Dbase = (
    802816*rho**4*w**4
    +476656*rho**2*w**6
    +41330688*rho**2*w**2
    +19393941*w**4
    +531948096
)/(14776336*w**2)

D = (
    Dbase
    +8*w**2*z**2/31
    +w**2*b**2/4
    +16*b*rho*w**2*z/31
)

completed = (
    Dbase
    +w**2*(b+32*rho*z/31)**2/4
    +8*w**2*(31-32*rho**2)*z**2/961
)

print("STAGE 1: square completion", flush=True)

assert S.simplify(D-completed) == 0
print("PASS: exact square completion")

print("STAGE 2: quadratic discriminant", flush=True)

poly = S.Poly(D, b)
disc = S.factor(S.discriminant(poly.as_expr(), b))

expected_disc = S.factor(
    -w**2*Dbase
    +8*w**4*(32*rho**2-31)*z**2/961
)

assert S.simplify(disc-expected_disc) == 0

print("PASS: exact discriminant")
print("DISCRIMINANT =", disc)

print("STAGE 3: recover R28U", flush=True)

special_b = -32*rho*z/31

special_D = S.factor(D.subs(b, special_b))

expected_special = (
    Dbase
    +8*w**2*(31-32*rho**2)*z**2/961
)

assert S.simplify(
    special_D-expected_special
) == 0

print("PASS: R28U specialization recovered")

print("STAGE 4: explicit quadratic branches", flush=True)

branches = [
    S.factor(-32*rho*z/31 + 2*S.sqrt(disc)/w**2),
    S.factor(-32*rho*z/31 - 2*S.sqrt(disc)/w**2)
]

for i, branch in enumerate(branches, 1):
    assert S.simplify(D.subs(b,branch)) == 0
    print("PASS: branch", i)

print("R29E PASS")
print("REAL DOMAIN: w != 0, rho^2 > 31/32")
print("REAL BRANCHES REQUIRE: discriminant >= 0")
print("GENERAL Ddot: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
