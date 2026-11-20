import sympy as S

print("=== CP098-T90C-R29J ===", flush=True)

w, rho, z, a, b = S.symbols(
    "w rho z a b", real=True
)

lam_general = (
    -w**2*(961*a+448*b*rho**2+448*rho*z)
    /(4*(112*rho**2*w**2+2883))
)

lam_reduced = (
    -112*w**2*rho*(b*rho+z)
    /(112*rho**2*w**2+2883)
)

print("STAGE 1: reduced multiplier")

assert S.simplify(
    lam_general.subs(a,0)-lam_reduced
) == 0

print("PASS: exact a=0 reduction")

print("STAGE 2: denominator domain")

den = 112*rho**2*w**2+2883

assert S.simplify(den.subs({w:1,rho:2})) > 0

print("PASS: denominator positive at witnesses")
print("For real w,rho: denominator >= 2883 > 0")

print("STAGE 3: R29G regressions")

wv = S.Integer(1)
rv = S.Integer(2)

Dbase = (
    802816*rv**4
    +476656*rv**2
    +41330688*rv**2
    +19393941
    +531948096
)/S.Integer(14776336)

z2min = S.factor(
    961*Dbase/(8*(32*rv**2-31))
)

zv = S.sqrt(z2min+1)

disc = S.Rational(776,961)

center = -S.Rational(32,31)*rv*zv

branches = {
    "PLUS": center+2*S.sqrt(disc),
    "MINUS": center-2*S.sqrt(disc)
}

expected = {
    "PLUS": (
        -S.Rational(1792,103261)*S.sqrt(194)
        +S.Rational(14,3201091)*S.sqrt(144209559530)
    ),
    "MINUS": (
        S.Rational(1792,103261)*S.sqrt(194)
        +S.Rational(14,3201091)*S.sqrt(144209559530)
    )
}

for name,bv in branches.items():
    actual = S.simplify(
        lam_reduced.subs({
            w:wv,rho:rv,z:zv,b:bv
        })
    )

    assert S.simplify(actual-expected[name]) == 0
    print("PASS:",name,"R29G multiplier regression")

print("R29J PASS")
print("MULTIPLIER FORMULA: REGRESSION VERIFIED")
print("GENERAL CANONICAL RE-DERIVATION: NOT PERFORMED")
print("Ddot FAMILY: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
