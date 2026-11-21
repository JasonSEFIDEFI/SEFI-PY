from pathlib import Path
import sympy as S

print("=== CP098-T90C-R29AS ===", flush=True)

folder = Path(__file__).parent
path = folder / "CP098_T90C_R29AQ.py"
source = path.read_text()

marker = '\nprint("STAGE 2: isolate coefficient roots"'
assert source.count(marker) == 1, "STOP: source boundary ambiguous"

ns = {
    "__name__": "r29as_definitions",
    "__file__": str(path),
}
exec(compile(source.split(marker)[0], str(path), "exec"), ns)

Px = ns["Px"]
u = ns["u"]

degrees = sorted(k for (k,), c in Px.terms() if c != 0)
lowest, highest = degrees[0], degrees[-1]

C0 = S.Poly(Px.nth(lowest), u)
CL = S.Poly(Px.nth(highest), u)

assert lowest == 0
assert highest == 26

assert all(c > 0 for c in C0.all_coeffs()), (
    "STOP: positive endpoint not certified"
)
assert all(c < 0 for c in CL.all_coeffs()), (
    "STOP: negative endpoint not certified"
)

print("PASS: exact positive constant coefficient")
print("PASS: exact negative leading coefficient")
print("CERTIFIED: Delta positive near x=0 for every u>0")
print("CERTIFIED: Delta negative for sufficiently large x")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("R29AS COMPLETE")
