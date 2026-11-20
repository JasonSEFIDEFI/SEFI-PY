
import sympy as S

print("=== CP098-T90C-R28U ===")

w, rho, z = S.symbols(
    "w rho z", real=True, nonzero=True
)

Dbase = (
    802816*rho**4*w**4
    +476656*rho**2*w**6
    +41330688*rho**2*w**2
    +19393941*w**4
    +531948096
)/(14776336*w**2)

b = -32*rho*z/31

D = (
    Dbase
    +8*w**2*z**2/31
    +w**2*b**2/4
    +16*b*rho*w**2*z/31
)

reduced = S.factor(D)

expected = (
    Dbase
    +8*w**2*(31-32*rho**2)*z**2/961
)

assert S.simplify(reduced-expected) == 0

z2 = S.factor(
    961*Dbase/(8*w**2*(32*rho**2-31))
)

assert S.simplify(
    reduced.subs(z**2,z2)
) == 0

print("PASS: exact reduced D identity")
print("PASS: D=0 algebraic solution")
print("DOMAIN: w != 0, rho^2 > 31/32")
print("z^2 =", z2)
print("R28U PASS")
print("FAMILY-LEVEL Ddot: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
