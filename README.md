# SEFI-PY

**Geometric modeling, classical field simulations, and reproducible research by Jason Duran Dutton.**

SEFI-PY is a research repository containing several independent Python and JavaScript projects. It includes geometric trajectory demonstrations, classical identity dynamics, a synthetic microwave diagnostic workbench, a local vortex research console, an educational double-slit explorer, and separate recovery benchmarks. It is not one installable Python distribution or a single validated physical theory.

SEFI (Single Entity Field Interpretation), DEFI (Dynamic Entity Field Integration), and GWFM are the author's research frameworks. Their terminology motivates some models and experiments; using that terminology does not establish their physical interpretation. Numerical consistency, engineering utility, and experimental validity are separate questions.

This guide describes the default-branch source at [`847a071`](https://github.com/JasonSEFIDEFI/SEFI-PY/tree/847a071351ed2b1c54e240fd9a520a1cce238178). Research on other branches and in open pull requests is not included in the capabilities below.

## Choose a starting point

| Goal | Project | What it actually does |
|---|---|---|
| Run a small geometric simulation | [Evolution engine](SEFI-PY%20Evolution/main.py) | Advances a 3D Frenet-frame trajectory, maps a DNA string to curvature, classifies geometric displacement, and renders voxel occupancy as ASCII. |
| Check classical identity mathematics | [Time dynamics](SEFI-PY/docs/TIME_DYNAMICS.md) | Computes prescribed-geometry identity equations, spectra, energy identities, directed clock phase, exact rational recurrence certificates, and Floquet diagnostics. |
| Reproduce measurement-to-diagnosis analysis | [Calibration workbench](SEFI-PY%20Evolution/research/calibration_workbench/README.md) | Simulates driven two-level populations and compares independent fitting, conventional EWMA tracking, and local tangent tracking; supports measurement-only CSV replay. |
| Inspect and evolve computed classical fields | [Vortex Research Console](SEFI-PY%20Evolution/sim/vortex_console/README.md) | Solves axisymmetric profiles, performs restricted physical evolution and separate fixed-charge relaxation, and stores runs and checkpoints. |
| Explore quantum interference visually | [Double Slit Explorer](SEFI-PY%20Evolution/sim/double_slit_explorer) | Samples a standard scalar diffraction/interference reference and presents interpretive chapters. |
| Evaluate auxiliary-probe recovery | [Displacement recovery benchmark](SEFI-PY/benchmarks/displacement_recovery/README.md) | Compares estimators under a specified shared-displacement pilot/payload model, with matched-information controls. |

For an engineering demonstration, start with the calibration workbench. For a small mathematical example with no numerical runtime dependency, start with Time dynamics.

## Repository architecture

These projects share research context, but do not form one execution pipeline.

```text
SEFI-PY/
  core/, sim/, utils/         Earlier Python identity/entity engine
  core/time_dynamics.py      Standalone classical identity mathematics
  tests/, demo_*.py           Python checks and demonstrations
  benchmarks/displacement_recovery/
                             Independent pilot-assisted recovery study
SEFI-PY Evolution/
  core/                      Four-coordinate state and diagonal metric
  geometry/                  Curves, transforms, Frenet-frame updates
  dna/                       Sequence encoding and curvature mapping
  sim/                       Geometric trajectory integration
  qec/                       Geometric displacement/projection helpers
  voxel/                     Discretization and ASCII projection
  research/calibration_workbench/
                             Independent NumPy diagnostic experiment
  research/physical_matter_*/
                             Dated classical-field research checkpoints
  sim/vortex_console/        Independent NumPy solver and loopback web UI
  sim/double_slit_explorer/   React/TypeScript/Vite educational app
research/                    Dated experiments, notes, and stored results
SEFI QEC Stack/              Separate QEC prototype, with nested copies
sefi-*/, sefi/, package.json Earlier JavaScript modules and CLI;
                             some directories also contain Python sketches
examples/                    Legacy integration examples
```

The Evolution demo calls sequence encoding → curvature mapping → trajectory integration → voxel rendering, with a separate illustrative displacement classification. Its `SEFISpace` metric is a diagonal quadratic form: `metric_distance` returns the squared interval, not its square root.

The workbench defines its reference in measurable population space; it does **not** import the geometric engine. The console likewise has its own equations and solver. The browser explorer uses a standard quantum reference rather than deriving interference from the Python models.

Run each Python project from its own directory. Several folders use generic imports such as `core`, `sim`, and `qec`; combining their import paths can select the wrong module. The root `package.json` runs the older Node CLI, not the React explorer or Python programs.

## Installation

Clone the repository and create a fresh environment. The repository contains historical environments and archives; do not rely on a checked-in `.venv` or ZIP as the current source.

```sh
git clone https://github.com/JasonSEFIDEFI/SEFI-PY.git
cd SEFI-PY
python -m venv .venv
```

Activate it with the command for your shell:

```powershell
# Windows PowerShell
.\.venv\Scripts\Activate.ps1
```

```sh
# macOS/Linux
source .venv/bin/activate
```

Use Python 3.12 as the previously documented reference environment for the console and benchmark. The workbench requires Python 3.10 or newer. The focused checks below were also rerun with Python 3.14.7, NumPy 2.5.2, and pytest 9.1.1; this is not a compatibility guarantee for every legacy module.

Install the dependencies for the project you choose. There is no root Python packaging configuration for the whole repository; `pip install sefipy` is not a documented installation route. Requirements are separate and generally not locked. Record Python and package versions with your results, and retain the exact source commit.

## Reproducible usage

All paths below start at the repository root unless stated otherwise. Use a fresh terminal at the root when switching projects.

### 1. Evolution engine

```sh
cd "SEFI-PY Evolution"
python -m pip install -r requirements.txt
python -m pytest tests -q
python main.py
```

The demo encodes `ATCGATCG`, generates 20 trajectory points, prints a geometric error classification, and displays an ASCII density projection. It injects unseeded random noise, so exact coordinates and artwork are not reproducible from the default command. For a controlled experiment, set and record the NumPy random seed before calling the solver, along with curvature, torsion, time step, and noise amplitude.

This is a geometric demonstration. In [`DEFIIntegrator.step`](SEFI-PY%20Evolution/sim/defi_integrator.py), correction uses `self.pos - noise`, giving it the actual injected perturbation reference. That is an oracle-assisted simulation, not a measurement-only decoder or evidence of quantum error correction. The Frenet update normalizes each vector separately; it does not explicitly enforce pairwise orthogonality after each step.

### 2. Time dynamics

The demonstration uses only the Python standard library; pytest is needed for its tests.

```sh
cd SEFI-PY
python -m pip install pytest
python demo_time_dynamics.py
python -m pytest tests/test_time_dynamics.py -q
```

The demo prints the exact frequency ratios `(1/5, 1, 7/5)`, recurrence after five Origin cycles (`10*pi`), a nonlinear barrier of `0.5`, and stable/unstable Floquet comparisons for the specified classical model. These are dimensionless model calculations. Geometry is prescribed; phase locking here is not an implemented attracting synchronization mechanism.

See [mathematical scope](SEFI-PY/docs/TIME_DYNAMICS.md) and the [recorded validation](SEFI-PY/docs/TIME_VALIDATION.md). The earlier Python engine has additional NumPy/plotting dependencies that its minimal `requirements.txt` does not fully describe; the commands above target only this standalone module.

### 3. Synthetic microwave calibration workbench

```sh
cd "SEFI-PY Evolution/research/calibration_workbench"
python -m pip install -r requirements.txt
python -W error -m unittest -v
python workbench.py --output demo_output --runs 40 --seed 20260925
python workbench.py --input demo_output/measurements.csv --output replay_output
```

Open `demo_output/index.html` locally. Outputs include benchmark metrics, synthetic measurements, and service reports; replay writes `replay_output/replay_report.md`. GitHub displays HTML source rather than hosting this local page.

The design uses 12 settings per scan, 1,024 shots per setting, fixed 2% symmetric readout error, and dimensionless frequencies/times relative to `Omega0`. Thresholds are calibrated on separate healthy simulations. Recovery requires three consecutive acceptable scans with residual checks. The program recommends checks and reports uncertainty; it does not apply a correction or control equipment.

CSV columns are `scan,setting,offset_omega0,duration_omega0,shots,excited`. Scans must start at zero, be contiguous and complete, and match the configured settings and shot count. Replay excludes truth labels but still uses the synthetic reference and calibration assumptions. See the [full schema and physical model](SEFI-PY%20Evolution/research/calibration_workbench/README.md).

### 4. Vortex Research Console

```sh
cd "SEFI-PY Evolution/sim/vortex_console"
python -m pip install -r requirements.txt
python -m unittest test_console test_rest_console -v
python -m unittest test_dispersion -v
python server.py
```

Open [http://127.0.0.1:8765](http://127.0.0.1:8765); stop the server with Ctrl+C. It binds to loopback and uses no external AI service. New runs are stored under this project's `runs/` directory. Retain parameters, model version, full checkpoints, and results when sharing a numerical claim.

The console implements coupled classical fields, Newton profile searches, and fourth-order Runge–Kutta evolution. Physical evolution is restricted to `lambda=25, N=0`. Fixed-charge relaxation is a separate computation with a separate checkpoint format. Replay and camera animation do not run new dynamics. Read the [equations and limits](SEFI-PY%20Evolution/sim/vortex_console/README.md), [fixed-charge research](SEFI-PY%20Evolution/sim/vortex_console/FIXED_CHARGE_RESEARCH.md), and [two-charge accounting](SEFI-PY%20Evolution/sim/vortex_console/TWO_CHARGE_ACCOUNTING.md).

### 5. Double Slit Explorer

Use Node.js 24 and npm, following the project's recorded installation guide and `.nvmrc`. This is a separate frontend installation.

```sh
cd "SEFI-PY Evolution/sim/double_slit_explorer"
npm install
npm test
npm run build
npm run dev
```

Open the local address printed by Vite. The app samples a normalized scalar Fraunhofer reference with visibility controlled by an idealized distinguishability model. All chapters use the same detector distribution; GWFM/SEFI/DEFI graphics are interpretations, not independent predictions. See [installation](SEFI-PY%20Evolution/sim/double_slit_explorer/INSTALL.md) and [validation](SEFI-PY%20Evolution/sim/double_slit_explorer/docs/VALIDATION.md).

### 6. Pilot-assisted displacement recovery

```sh
cd "SEFI-PY/benchmarks/displacement_recovery"
python -m pip install -r requirements.txt
python -m unittest -v test_benchmark
python run_benchmark.py
python make_figures.py
python build_report.py
```

Use Python 3.12 or newer. Simulation uses NumPy, plots use Matplotlib, and PDF generation uses ReportLab. Read the [protocol and reproduction guide](SEFI-PY/benchmarks/displacement_recovery/README.md) before interpreting the results. Sixteen auxiliary coherent probes share the payload displacement; this is auxiliary-noise sensing, not a stabilizer-code threshold experiment. The geometry estimator coincides with a conventional Gaussian posterior mean. An improvement over maximum likelihood under that prior is not a SEFI-specific or quantum advantage.

## Validation status

The following distinguishes checks rerun for this README from historical reports. It is not a repository-wide certification.

| Component | Evidence | Status and scope |
|---|---|---|
| Evolution engine | Four tests in its `tests/`; `main.py` | **Rerun for this README:** 4 passed; demo completed. These are small metric, geometry, displacement, and integrator checks. |
| Time dynamics | `tests/test_time_dynamics.py`; standalone demo | **Rerun:** 12 passed; demo completed. Tests compare analytic and integrated model behavior, energy identities, recurrence, phase sampling, and Floquet calculations. |
| Calibration workbench | `test_workbench.py`; seeded benchmark and CSV replay | **Rerun:** 11 passed; 40 runs per scenario and replay completed. Tests include physical limits, observability, malformed input, recovery gates, and repeatability. |
| Earlier Python engine, full suite | [Time validation record](SEFI-PY/docs/TIME_VALIDATION.md) | Historical report: 33 tests passed under Python 3.12. The full suite was not rerun for this README. |
| Vortex console | [Test report](SEFI-PY%20Evolution/sim/vortex_console/TEST_REPORT.md) | Historical numerical/API and limited browser checks; includes checkpoint continuation and refinement. Not rerun here. Charge-tab visual verification remained outstanding in that report. |
| Double-slit explorer | [Validation record](SEFI-PY%20Evolution/sim/double_slit_explorer/docs/VALIDATION.md) and source tests | Historical production build and 7 numerical tests passed. Browser coverage is partial; mobile, reduced-motion, full chapter coverage, and image download checks were not all completed. Not rerun here. |
| Recovery and dated research studies | Project-specific protocols, tests, reports, and stored outputs | Inspect each experiment's assumptions and revision. Not rerun for this README. |

The reruns used isolated copies of source from `847a071`, Windows, Python 3.14.7, NumPy 2.5.2, and pytest 9.1.1. No browser or hardware validation was performed for this documentation update.

For the seeded workbench rerun, healthy per-scan false alarms were **1.11%** (independent fit), **0.94%** (conventional EWMA), and **1.22%** (tangent EWMA). Each method flagged all 120 injected-fault runs at least once within the fault windows. Median detection delays matched across methods; smoothing lagged during recovery. These reproduce the [recorded benchmark](SEFI-PY%20Evolution/research/calibration_workbench/VALIDATION.md). **No advantage over conventional tracking is established**, and these synthetic frequencies are not hardware reliability guarantees.

## Scientific and engineering limitations

- **Classical model results are conditional.** Finite-grid solutions, small residuals, conserved model charges, and passing tests do not establish continuum existence, full nonlinear stability, observed particles, quantum statistics, universal matter coupling, or Einstein gravity.
- **Vortex evidence is restricted.** The main evolution is finite-domain and axisymmetric. The documented rest-state obstruction applies to its specified model window. A different parameter window or a charged-torus relaxation candidate does not prove stable matter; the candidate lacks an original-field vortex zero. Its constrained energy is not a calibrated particle mass.
- **Synthetic diagnostics omit real-device effects.** Same-model generation and inference favor the estimators. Leakage, multilevel dynamics, colored noise, instrument latency, uncertain references, and correlated drift need additional models and data. Local uncertainty estimates are not validated confidence intervals.
- **Research labels are not hardware evidence.** DNA-to-curvature mappings do not establish a biological mechanism. Geometric displacement flags are not measured quantum syndromes. The serial interface is a sketch, not an experimentally validated transducer system.
- **Legacy examples need integration work.** Root `examples/demo_sefi_worldline.py` imports underscore package names while the corresponding root directories use hyphens, supplies random placeholder transducer geometry, and attempts serial output to `COM3`. It is not a plug-and-play quickstart.
- **No AI efficiency breakthrough is demonstrated.** The [identity-memory exploration](research/2026-10-18/sefi-ai-exploration/FINDINGS.md) found no advantage over ordinary indexed state storage. It does not implement a replacement language model, establish token elimination, or measure data-center energy savings.
- **No propulsion result is established.** The [candidate field-recoil simulator](research/2026-10-18/sefi-dynamics/README.md) computes an explicitly hypothetical 1D interaction. Its browser controls play stored results; source-to-material coupling and experimental calibration remain missing.

## Research context and reading order

Evaluate implemented mathematics and measurements without assuming the wider interpretation is correct.

1. Read the selected project's model, units, protocol, and tests.
2. Read its negative controls, failed searches, and validation limits alongside any favorable result.
3. Retain the exact commit and distinguish stored output from a newly executed calculation.
4. Consult the foundational drafts separately: [SEFI-Manuscripts](https://github.com/JasonSEFIDEFI/SEFI-Manuscripts) and the [October research index](RESEARCH-UPDATE-2026-10-18.md). Drafts and submission history are not evidence of journal acceptance or experimental confirmation.

For field research, start with the [23 September checkpoint](SEFI-PY%20Evolution/research/physical_matter_2026_09_23/README.md), the [console reviewer guide](SEFI-PY%20Evolution/sim/vortex_console/REVIEWER_GUIDE.md), and the [19 October research directory](research/2026-10-19). Later evidence includes interrupted runs and unresolved angular convergence; the console's evidence page must not be read as proof of general 3D stability.

The workbench is an independent engineering proposal using synthetic data, not an eleQtron assignment or endorsed product. Related [SEFI_QEC_Stack](https://github.com/JasonSEFIDEFI/SEFI_QEC_Stack) work is a separate repository; its capabilities and validation do not automatically transfer here.

## Roadmap — proposed work, not implemented guarantees

| Direction | Next discriminating milestone |
|---|---|
| Packaging and reproducibility | Separate project environments, complete dependency manifests, locked reproduction environments, and clean-checkout CI for explicitly supported entry points. |
| Measurement-based diagnostics | Learn a reference from approved hardware data, test mismatch and drift outside the synthetic generator, and compare against existing methods at equal information and measurement cost. |
| Recovery evaluation | Expand bias/mismatch controls and report when prior-assisted methods fail; retain a matched conventional posterior baseline. |
| Field mathematics | Resolve charge accounting and relevant existence/stability questions; broaden domain, grid, time-step, and genuinely 3D perturbation studies before making physical identifications. |
| UI verification | Complete documented browser, accessibility, mobile, export, and Charge-tab checks. |
| Identity memory | Compare ordinary entity memory with and without proposed identity constraints at equal answer quality, evidence retention, and total processing cost. |

These are directions for further work, not a release schedule. Research branches and open pull requests require their own review before they can be described as default-branch features.

## Contributing and citing

Report the project path, source commit, Python/Node and package versions, exact command, seed, parameters, expected behavior, and observed output. For numerical claims, include raw results and the controls needed to distinguish a model effect from discretization or estimator assumptions. Keep unrelated projects and their checkpoints separate.

Cite Jason Duran Dutton, this repository URL, the exact commit, the specific experiment or module, and the relevant manuscript when using substantial research material. See the [research data policy](research/RESEARCH_DATA_POLICY.md) for provenance and handling of research artifacts.

## License

Read the root [SEFI-PY Research License v1.0](LICENSE): it permits attributed research, educational, scientific, and other non-commercial use; commercial use requires prior written permission. Derivative works must identify modifications.

Some older package metadata declares MIT, and bundled subprojects contain their own license files. That inconsistency is not resolved by this README. Check the applicable files and obtain clarification from the author before relying on conflicting metadata.
