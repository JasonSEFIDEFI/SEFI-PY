
import sympy as S
import runpy, contextlib, io
from pathlib import Path

print("=== CP098-T90C-R28Q ===", flush=True)

source = Path(__file__).with_name("CP098_T90C_R28K2.py")
with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

base = ns["ns"]
point = ns["point"]
at = ns["at"]
bracket_at = ns["bracket_at"]

w = S.Integer(1)
rho = S.Integer(2)

Dbase = (
    802816*rho**4*w**4
    +476656*rho**2*w**6
    +41330688*rho**2*w**2
    +19393941*w**4
    +531948096
)/(14776336*w**2)

z2 = S.factor(
    Dbase*961/(8*w**2*(32*rho**2-31))
)
assert z2 > 0
z = S.sqrt(z2)
b = -32*rho*z/31

IS0 = 12/w**2 + S.Rational(448,961)*rho**2
M0 = 1+(S.Rational(9,4)+IS0**2+1+b**2)/2

point[base["q"][1]] = w
point[base["P"][2]] = rho
point[base["p1"][2]] = z
point[base["p1"][3]] = 0
point[base["IS"]] = IS0
point[base["pA"]] = 0
point[base["pS"]] = 0
point[base["pW"]] = b
point[base["p0"][0]] = 2*w**2*rho**2/31-M0
point[base["p0"][1]] = -8*w*rho*z/31

print("STAGE 1: family formula", flush=True)
Dfamily = (
    Dbase
    +8*w**2*z**2/31
    +w**2*b**2/4
    +16*b*rho*w**2*z/31
)
assert S.simplify(Dfamily) == 0
print("PASS: exact family D=0", flush=True)

print("STAGE 2: canonical verification", flush=True)
for name in (
    "phi1", "phi2", "chi", "psi", "theta",
    "B", "A", "J", "Q", "pA"
):
    value = at(base[name])
    print(name, "=", value, flush=True)
    assert value == 0

Dcanonical = bracket_at(base["B"], base["H"])
print("D canonical =", Dcanonical, flush=True)
assert Dcanonical == 0

print("R28Q PASS: SECOND EXACT D=0 WITNESS")
print("HIGHER TANGENCY: NOT YET TESTED")
print("FULL DIRAC CLOSURE: OPEN")
