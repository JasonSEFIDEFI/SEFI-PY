# Vortex Research Console

An isolated local research instrument for inspecting computed field profiles, testing axisymmetric field evolution, and exploring fixed-charge relaxation. **Research prototype, console version 0.2.** The existing SEFI-PY engine and SEFI-QEC are unchanged. This is a locally run application, not a hosted online simulation.

## New in 0.2: fixed-charge research

The **Charge** tab loads the new charged-torus candidate and failed lower-charge searches, or executes an energy-decreasing fixed-carrier-charge search. It displays relaxation iteration, constrained energy per charge, carrier frequency, original-field minimum amplitude and residual. It supports pause/resume/stop, saved frames, checkpoint download and deterministic checkpoint continuation. Its checkpoint format is separate from physical-evolution checkpoints.

Run `python -m unittest test_console test_rest_console -v` to check both computational modes. The moving-ring equations and model version remain unchanged. The added `rest_search.py` supplies the research computation without modifying the established engine. Read [FIXED_CHARGE_RESEARCH.md](FIXED_CHARGE_RESEARCH.md) and [TWO_CHARGE_ACCOUNTING.md](TWO_CHARGE_ACCOUNTING.md) before interpreting the new branch.

For independent command-line reproduction, run `python rest_search.py --Q 1000 --lam 12 --N 1 --L 24 --h .5 --steps 18000`. Use `--initial Q1000_wide --h .4` to start a refinement from the saved wider-domain candidate. Results are written into a new unique folder under `runs/`. Run `python two_charge_audit.py` to reproduce the compact-window charge-compensation calculation; its result is written to `runs/two_charge_torus_audit.json`.

The Q=1000 candidate has a toroidal carrier but no original vortex zero. Its convergence does not establish stability. The Charge mode holds carrier charge fixed and original background frequency fixed; it does not conserve original-field excess charge during relaxation. Compensation tests account for both charges only in prepared initial data. They do not constitute a new physical evolution. The displayed energy is not particle mass. Current physical evolution remains restricted to the original lambda=25,N=0 moving-ring model.

## Open the console

With Python 3.12 and NumPy installed, run `python server.py` from this folder, then open **http://127.0.0.1:8765**. Windows users can run `launch.ps1`; it uses the available bundled Python runtime when present, otherwise the Python on PATH. Keep the server running while using the console. Ctrl+C stops it. The console binds only to loopback and makes no external AI calls. Rendering needs no downloaded JavaScript libraries.

Tested here with Python 3.12.14 and NumPy 2.3.5. `requirements.txt` describes the dependency range; other versions in that range have not all been tested. Run `python test_console.py` for the numerical and local API checks.

## Explore

- **View:** orbit with a drag, zoom with the wheel, choose top/side views, show either field, change display thresholds and opacity, compare the saved baseline, and inspect density or phase in a meridional cross-section. Export a scene PNG or view-settings JSON.
- **Physics:** edit model coefficients, domain, grid, and seed radius; select **Calculate profile**. A Newton solver recomputes both fields. A seed radius is an initial guess, not a prescribed result. Changing a parameter alone does not morph the current display.
- **Evolve:** select a disturbance and time settings, then **Execute evolution**. Pause/resume/stop, save a checkpoint, replay stored frames, or continue a checkpoint to a later end time. Sweeps run independent fixed-parameter simulations for up to five disturbance, travel-coefficient, or carrier-frequency values.
- **Math:** see executable options, equations, geometric interpretations, assumptions, and proposed directions that are not implemented as dynamics.
- **Review:** copy questions tailored to general, theoretical, mathematical, experimental, numerical, or software review. See [REVIEWER_GUIDE.md](REVIEWER_GUIDE.md) for the full question set and Copilot setup.

**Replay** displays computed saved frames; it does not start a simulation. **Auto-orbit** moves the camera only. PNG output is a visualization, not a numerical checkpoint. Exported view settings include pending and displayed parameters separately.

## Mathematics implemented

Use signature (+---), model units with underlying Minkowski propagation speed 1, and complex fields with kinetic normalization 1/2. The potential is

`V = |Φ|²/2 + |Φ|⁴/4 + λ|Σ|⁴/4 + (4|Φ|²−3)|Σ|²/2`.

For the axisymmetric moving profile, `Φ=exp(i√2 t)ψ(r,Z)` and `Σ=exp(iνγ(t−vz))s(r,Z)exp(iNφ)`, where `Z=γ(z−vt)`, `γ=(1−v²)^−1/2`, and `v=c/√(8+c²)`:

`Δψ + ic ∂Zψ + (1−|ψ|²−4s²)ψ = 0`

`Δs − N²s/r² + (3+ν²−4|ψ|²−λs²)s = 0`.

Here Δ is the axisymmetric cylindrical Laplacian. The outer boundary is ψ=1, s=0; regularity holds at r=0, with s=0 there for N>0. Profile searches assume real ψ and s even across Z=0 and imaginary ψ odd. Full time evolution does not enforce this Z reflection symmetry after initialization.

`physics.py` implements a coupled Newton solve and fourth-order Runge–Kutta evolution. Evolution currently supports **λ=25, N=0 only**. With `T=t`, the moving-coordinate Laplacian uses spatial coordinates r,Z; the equations include mixed time/space advection and factored-phase terms. They are listed in the Math menu and in `Evolution.rhs`.

Carrier charge in evolution is integrated with the finite-volume radial weight:

`Q = ∫[ν|s|² + Im(s* ∂T s)/γ − v Im(s* ∂Z s)] 2πr dr dZ`.

The stored profile report uses a different quadrature: baseline Q=7.1333803996 versus the evolution control Q=7.1334944088. This small discretization difference is not a time-evolution drift. Adding the initial perturbation changes the initial charge; each evolution measures drift relative to its own initial Q.

Phase-aligned carrier distance minimizes the global phase offset using the weighted overlap with the reference carrier. Core winding is measured on a meridional contour around the initial ring radius; it is not carrier winding N. An undefined or invalid contour is reported without a winding value.

## Rendering and limits

The 3D surfaces are obtained by linearly interpolating density contours on the numerical r,Z grid and revolving them around the symmetry axis. There is no analytic torus substituted for calculated data. Surface opacity and thresholds affect visibility only. Rendering uses depth-sorted canvas polygons; overlapping transparent surfaces are illustrative and can show sorting artifacts. Cross-sections supply a direct alternative view of the sampled data. Density colors can clip above the displayed scale. Phase at an exact field zero is undefined.

The default scene shows the inner ±10 units in cross-section; the solver domain is larger. Original-field surface colors encode phase approximately; the cross-section phase option supplies the full cyclic phase view. Carrier azimuthal phase for N>0 is not rendered. Coordinates and times are dimensionless; no calibration to SI particle scales is established.

A low residual establishes a discrete solution only. Evolution is finite-domain and axisymmetric, with fixed outer boundaries; it excludes general 3D disturbances and may develop boundary reflections. The original wider-domain study showed roughly 3.9% ring-radius sensitivity. Short-time conservation and refinement tests do not prove nonlinear stability. At λ=25, the documented positivity/scaling argument obstructs the specified nontrivial time-harmonic rest profile. Removing that obstruction in another parameter window does not prove existence.

**No identification with observed physical matter, quantum statistics, universal coupling, or Einstein gravity has been derived.** Effective phase geometry is an approximation discussed in the menu, not a spacetime simulation.

## Reproducibility and future versions

Baseline data are retained in `data/`. New runs store parameters, model version, solved profiles, display frames, results and full checkpoints under `runs/<id>/`. Checkpoints include both complex fields, both velocities, time, step count, initial charge normalization, and the reference profile. Display JSON is rounded to six decimal places; checkpoints preserve numerical precision. Server restarts retain checkpoint files; the current-session job list is not a persistent history browser.

Runtime excludes pauses; stop/cancel during a profile solve takes effect between Newton iterations. Limits are bounds on local resources, not guarantees of numerical accuracy. Frame storage reserves space for the checkpoint; actual disk usage also includes profile and metadata files. Archive run folders if total storage exceeds the server limit. Increasing grid size can make dense Newton solves substantially more expensive.

Read [MODEL_HISTORY.md](MODEL_HISTORY.md) for the user's requested continuing-development policy: preserve the baseline, version equations and saved states, add discriminating tests, and update explanations as new mathematics is established. The console rejects checkpoints from an incompatible model version.

[TEST_REPORT.md](TEST_REPORT.md) records executed checks. [REVIEWER_GUIDE.md](REVIEWER_GUIDE.md) contains questions for reviewers and GitHub Copilot usage. Repository-wide Copilot guidance lives at the SEFI-PY repository root in `.github/copilot-instructions.md`. Guidance does not share the author's personal Copilot session or guarantee an accurate answer.
