import sympy as S

print("=== CP098-T90C-R29F ===", flush=True)

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

disc = S.factor(S.discriminant(D,b))

z2min = S.factor(
    961*Dbase/(8*w**2*(32*rho**2-31))
)

assert S.simplify(
    disc.subs(z**2,z2min)
) == 0

print("PASS: exact discriminant threshold")

print("STAGE 2: construct off-center witness")

wv = S.Integer(1)
rv = S.Integer(2)

threshold = S.factor(
    z2min.subs({w:wv,rho:rv})
)

zv = S.sqrt(threshold+1)

disc_v = S.simplify(
    disc.subs({w:wv,rho:rv,z:zv})
)

assert disc_v > 0

b_center = -S.Rational(32,31)*rv*zv

b_plus = S.simplify(
    b_center + 2*S.sqrt(disc_v)/wv**2
)

b_minus = S.simplify(
    b_center - 2*S.sqrt(disc_v)/wv**2
)

assert S.simplify(b_plus-b_minus) != 0

for branch in (b_plus,b_minus):
    assert S.simplify(
        D.subs({
            w:wv,rho:rv,z:zv,b:branch
        })
    ) == 0

print("PASS: two distinct real D=0 branches")
print("w =",wv)
print("rho =",rv)
print("z^2 =",S.simplify(zv**2))
print("discriminant =",disc_v)
print("b_plus =",b_plus)
print("b_minus =",b_minus)

print("R29F PASS")
print("OFF-CENTER D=0 WITNESS: CONSTRUCTED")
print("Q TANGENCY: NOT YET RECHECKED")
print("Ddot: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
