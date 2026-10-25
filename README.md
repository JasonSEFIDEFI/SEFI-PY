# SEFI-PY — field engineering and reproducible research

## Double Slit Explorer — educational demo

[Open the Double Slit Explorer project](SEFI-PY%20Evolution/sim/double_slit_explorer) for a five-chapter introduction to quantum interference, QFT, and the proposed GWFM, SEFI and DEFI interpretations. Includes interactive detection patterns, observation controls, identity cards and a geometry–identity illustration. Run `npm install` then `npm run dev` inside that folder. [Installation guide](SEFI-PY%20Evolution/sim/double_slit_explorer/INSTALL.md) · [Validation](SEFI-PY%20Evolution/sim/double_slit_explorer/docs/VALIDATION.md).

This is an educational app separate from the Vortex Research Console. All detector probabilities use the same standard quantum reference; the proposed geometry does not independently derive matter or interference. Production build and seven numerical tests passed; browser coverage is partial and documented.

Python research software by **Jason Duran Dutton**, connecting geometric modeling with practical questions about measurement, drift, diagnosis, and recovery.

## For the meeting with Hillary — Friday, 25 September 2026

The proposed contribution is a small, testable engineering workflow: turn calibration measurements into evidence of a change, identify the next useful check, and document whether the system has recovered. Jason brings field experience in electrical integration, controls, commissioning, maintenance, and fault isolation; this prototype demonstrates how that approach can be extended through software.

### Start with the microwave calibration workbench

[**Open the prototype and reproduction guide**](SEFI-PY%20Evolution/research/calibration_workbench/README.md)

| What to review | Evidence |
|---|---|
| What was built and how to run it | [Workbench README](SEFI-PY%20Evolution/research/calibration_workbench/README.md) |
| What an engineer would receive | [Sample service report](SEFI-PY%20Evolution/research/calibration_workbench/sample_output/service_report.md) |
| When the tool declines to infer a cause | [Ambiguous measurement example](SEFI-PY%20Evolution/research/calibration_workbench/sample_output/ambiguous_report.md) |
| What the benchmark found | [Validation record](SEFI-PY%20Evolution/research/calibration_workbench/VALIDATION.md) and [metrics](SEFI-PY%20Evolution/research/calibration_workbench/sample_output/metrics.csv) |
| Implementation and checks | [Python source](SEFI-PY%20Evolution/research/calibration_workbench/workbench.py) and [tests](SEFI-PY%20Evolution/research/calibration_workbench/test_workbench.py) |

For the visual demonstration, download the workbench folder and open `sample_output/index.html` locally. GitHub does not render that HTML as a hosted application.

### The engineering question

Can a reproducible diagnostic workflow distinguish microwave amplitude drift from detuning, recognize an inadequate measurement design, and verify restoration against explicit acceptance criteria?

The motivation comes from [eleQtron researchers' public work on microwave-chain errors and coherent storage transfer](https://indico.fysik.su.se/event/9371/contributions/15175/) and the company's [RF/system-test responsibilities](https://eleqtron.com/en/jobs/rf-radio-frequency-test-engineer-f-m-d/). This is an independent proposal; it is not an eleQtron assignment or statement about an internal fault.

### What has been demonstrated

- A standard two-level quantum response generates synthetic calibration counts.
- Conventional fitting, conventional temporal tracking, and SEFI-inspired tangent tracking receive identical measurements, without simulation truth labels.
- Separate healthy calibration runs set alert thresholds. Held-out tests cover healthy operation and three drift scenarios.
- The tool reports ambiguous evidence and requires repeated acceptable measurements for recovery verification.
- Eleven tests pass for this addition, covering physical limits, input validation, observability, reproducibility, and recovery checks.

**The benchmark does not establish an advantage over conventional tracking.** Observed healthy per-scan false alarms are 1.11% for independent fitting, 0.94% for conventional tracking, and 1.22% for the SEFI-inspired method. All three detect all 120 simulated fault runs at least once during the injected fault windows. Smoothed methods lag during recovery. These are small synthetic experiments, not hardware reliability guarantees.

### Where the research framework enters

SEFI/DEFI motivates tracking deviations from a reference geometry. This prototype defines that geometry through a physical response model and measurable populations. Its weighted tangent projection and EWMA are established mathematical techniques; their inclusion does not establish a novel algorithm or validate the wider field theory. It neither implements quantum error correction nor controls equipment.

### Proposed next step with eleQtron

Identify one useful calibration workflow with an RF validation or operations engineer. Agree on the measurement format, relevant operating limits, and a success criterion. Evaluate approved anonymized data against the team's current method, keeping whichever method offers the strongest evidence. Potential measures include false alarms, time to diagnosis, recovery verification, and measurement burden.

## Broader research portfolio

The engineering prototype is separate from the foundational research below. Readers can evaluate its code and measurements without accepting the wider theoretical interpretation.

- [Vortex Research Console](SEFI-PY%20Evolution/sim/vortex_console/README.md): interactive views, profile recalculation, evolution, parameter sweeps, and saved checkpoints.
- [Console reviewer guide](SEFI-PY%20Evolution/sim/vortex_console/REVIEWER_GUIDE.md), [test report](SEFI-PY%20Evolution/sim/vortex_console/TEST_REPORT.md), and [model history](SEFI-PY%20Evolution/sim/vortex_console/MODEL_HISTORY.md).
- [23 September field research checkpoint](SEFI-PY%20Evolution/research/physical_matter_2026_09_23/README.md): coupled classical vortex/carrier calculations, limited perturbation tests, a rest-state obstruction, and unsuccessful revised searches.
- [Earlier evolution engine](SEFI-PY%20Evolution/README.md) and [QEC research repository](https://github.com/JasonSEFIDEFI/SEFI_QEC_Stack).
- [Professional portfolio and resume](https://github.com/JasonSEFIDEFI/CV-) and [doctoral research portfolio](https://github.com/JasonSEFIDEFI/PhD).

The foundational work investigates SEFI, DEFI, and GWFM proposals. Finite-domain classical solutions and passing software checks do not establish observed matter, a continuum existence theorem, nonlinear stability, universal gravity, or hardware QEC performance. Validation statements apply only to the named experiment and revision; this README does not certify all legacy modules.

## Research synchronization - October 18, 2026

- [Manuscript archive and author-review drafts](https://github.com/JasonSEFIDEFI/SEFI-Manuscripts/blob/main/RESEARCH-UPDATE-2026-10-18.md).
- [New simulations, assumptions, and reproducible calculations](https://github.com/JasonSEFIDEFI/SEFI-PY/tree/main/research/2026-10-18).
- [Book manuscripts and early source material](https://github.com/JasonSEFIDEFI/Publishing-V-Drive-Content).

These are research drafts and computational studies, not confirmations of new physics, a working propulsion engine, journal acceptance, or a doctoral award. Prior versions and submission history remain distinct. Existing draft pull requests are retained for review.

## Localized vacuum branch — 25 October 2026

[Continuation, binding threshold and publication-readiness record](research/2026-10-25/localized-branch/README.md): three-grid charge continuation, selected fragmentation comparisons, independent continuum collocation and an exploratory binding crossing near Q=1225.43 in the dimensionless lambda=12 model. Includes unsuccessful runs, raw profiles and limitations. **Exploratory / Not Promoted**; no observed-particle identification or general nonlinear stability claim.

## Stationary charge turn and critical coupling — 25 October 2026

[Frequency continuation and critical-limit study](research/2026-10-25/localized-frequency/README.md) separates a numerical charge minimum near Q=953.31 from the binding crossing near Q=1225.43 at lambda=12, compares two Q=1000 stationary states with different radial constrained-energy character, and tests a formal critical scaling problem as lambda approaches 16. Includes the stationary virial obstruction's assumptions, raw profiles, rejected continuation attempts and a publication assessment. **Exploratory / Not Promoted**; no general nonlinear stability or established novelty claim.

## First correction and working manuscript — 25 October 2026

[Critical-limit correction and author-review draft](research/2026-10-25/critical-correction/README.md) derives first-order charge coefficients without fitting full-model charge data, checks two declared closer couplings, and supplies a conditional amplitude-elimination bound. Prediction errors improve by factors of roughly 1,200–2,400 in those tests. Includes the [working manuscript](research/2026-10-25/critical-correction/MANUSCRIPT_DRAFT.md), a prior-work comparison, raw profiles and all numerical checks. **Exploratory / Not Promoted**; no uniform coupled-limit theorem or established novelty claim, and no journal submission.
