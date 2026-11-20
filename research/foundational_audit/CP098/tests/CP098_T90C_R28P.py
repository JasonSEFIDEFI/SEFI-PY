import sympy as S

print("=== CP098-T90C-R28P ===")

r = S.sqrt(S.Integer(1187904394))
F = -S.Rational(78415527239,7100029448)*r
K1 = S.Integer(0)
K2 = S.Rational(74324148,29791)
lam2 = S.Rational(7,2878195)*r

residual = S.simplify(F + lam2*K2)

expected = (
    -S.Rational(234725477359877,21264588196760)*r
)

assert residual == expected
assert residual != 0
assert K1 == 0

print("PASS: independent exact arithmetic")
print("Ddot =", residual)
print("PASS: Ddot is nonzero")
print("PASS: lambda1 cannot cancel obstruction")
print("R28J WITNESS: EXCLUDED BY HIGHER TANGENCY")
print("GENERAL D=0 BRANCH: STILL OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
