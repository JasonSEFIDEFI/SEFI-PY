import sympy as S

print("=== CP098-T90C-R29C ===", flush=True)

x, rho, z = S.symbols("x rho z", real=True)

z2 = (
    802816*rho**4
    + 41807344*rho**2
    + 551342037
) / (123008*(32*rho**2-31))

F = -2*rho*z*(
    154140672*rho**6
    + 9304325120*rho**4
    + 139108640128*rho**2
    - 70151578681
) / (28629151*(32*rho**2-31))

K2 = 12*(
    13488*rho**4
    + 639065*rho**2
    + 5541126
) / (29791*(32*rho**2-31))

b = -S.Rational(32,31)*rho*z

lam2 = -(
    448*b*rho**2 + 448*rho*z
) / (4*(112*rho**2+2883))

P = (
    17263755264*x**4
    + 1477480546304*x**3
    + 41987231391232*x**2
    + 389912088162608*x
    - 198668586495531
)

residual = S.factor(F + lam2*K2)

expected = -2*rho*z*P.subs(x,rho**2)/(
    28629151
    *(32*rho**2-31)
    *(112*rho**2+2883)
)

assert S.simplify(residual-expected) == 0
print("PASS: exact Ddot factorization", flush=True)

boundary = S.Rational(31,32)

assert S.simplify(
    P.subs(x,boundary)
) == 219821236558923

assert all(
    S.expand(P).coeff(x,k)>0
    for k in range(1,5)
)

print("PASS: P strictly increasing for x >= 31/32")
print("PASS: P(31/32) =", P.subs(x,boundary))

rho_w = S.Integer(2)
z_w = S.sqrt(S.factor(z2.subs(rho,rho_w)))

observed = S.simplify(
    residual.subs({rho:rho_w,z:z_w})
)

previous = (
    -S.Rational(
        2131753744725653,
        111262356920464396
    ) * S.sqrt(S.Integer(141894794986))
)

assert S.simplify(observed-previous) == 0
print("PASS: exact R28S/R28T regression")

print("R29C PASS")
print("RESULT: Ddot != 0 on specified D=0 subfamily")
print("SIX-PARAMETER FAMILY: NOT EXCLUDED")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
