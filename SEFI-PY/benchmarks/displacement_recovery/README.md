# Pilot-assisted bosonic displacement recovery benchmark

Jason D. Dutton | Preliminary computational evidence for NSF QIS | October 6, 2026

One experiment, three figures, complete raw/aggregate data and matched-information controls. This is a benchmark report, not a journal manuscript. It tests quantum-state recovery with a geometry-informed estimator motivated by reference distances. That estimator is mathematically a conventional Gaussian posterior mean. It does not validate SEFI, DEFI, GWFM, or a quantum advantage.

Start with `Benchmark_Report.pdf` or the editable `BENCHMARK_REPORT.md`. Read `TECHNICAL_AUDIT.md` for what already exists in the pinned SEFI-PY repository. Primary results are in `results/table_1_all_results.csv`; all three PNG/SVG figures are in `figures/`.

## Reproduce

Use Python 3.12 or later in a clean virtual environment. Install the versions in requirements.txt, then run:

```text
python -m pip install -r requirements.txt
python -m unittest -v test_benchmark
python run_benchmark.py
python make_figures.py
python build_report.py
```

The simulation requires NumPy only; plotting requires Matplotlib and the report requires ReportLab. ReportLab includes a binary wheel; no TeX installation is needed. The independent report-generating code is included. Re-run the repository audit from the original pinned checkout only; the audit script requires the parent repository and its inventory. It is not needed to reproduce the simulation.

Protocol: 20 seeds, 4,000 held-out trials per seed at each of seven noise levels; 560,000 total held-out trials, plus 280,000 calibration injections. Decoder A is exact maximum likelihood for the probe measurement likelihood; Decoder B is regularized phase-space geometry. A conventional prior-matched posterior is mandatory. Both principal decoders receive identical probe observations and have access to the same calibration information. The quantum payload is never measured to produce correction estimates.

The pilot-assisted assumption matters: sixteen coherent probe modes share a displacement with the unknown payload. This is an auxiliary-noise-sensing recovery study, not a stabilizer QEC code or threshold demonstration. Leakage is counted as failure; no postselection or normalized-projection fidelity is used.

The conventional posterior and geometry estimator coincide to numerical precision. An advantage over ML is therefore ordinary prior-assisted estimation, not evidence of a special geometric mechanism. Unmodeled displacement bias sharply degrades the prior-based decoder; see the report and stress results.

## Files

```text
displacement_recovery/
  protocol.json                 fixed design and primary endpoint
  bosonic.py                    oscillator reference and CPTP checks
  decoders.py                   ML, geometry, conventional and controls
  statistics.py                 seed-cluster bootstrap
  run_benchmark.py              simulation, controls, timing, quadrature
  test_benchmark.py              physical and inference tests
  audit_repository.py            static repository audit
  repository_inventory.json     pinned owned-tree hashes
  make_figures.py                Figures 1-3
  build_report.py                PDF from report source
  requirements.txt              exact installed package versions
  README.md
  TECHNICAL_AUDIT.md
  BENCHMARK_REPORT.md
  REPRODUCIBILITY.md
  NEGATIVE_CONTROLS.md
  IMPLEMENTATION_AND_PR.md
  NSF_PRELIMINARY_WORK.md
  Benchmark_Report.pdf
  audit/                        full source inventory and diagnostics
  results/                      raw NPZ, seed metrics, CIs, timing, checks
  figures/                      three PNGs and three SVGs
```

Use the existing repository research license and attribution. No changes to legacy modules are required. Suggested commits and a review-ready PR description are in `IMPLEMENTATION_AND_PR.md`; this package does not itself publish a GitHub branch or PR.
