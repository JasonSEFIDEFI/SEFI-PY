# Resumed research checkpoint — 2026-10-19

This checkpoint preserves completed follow-up calculations separately from the original interrupted run. All values are dimensionless. No main-console physical equations have been changed and no independent Authorship coupling has been introduced.

## Nonlinear angular checks

The new 24-angle rerun reached T=12 with dt=0.025. The 16-angle temporal rerun reached T=12 with dt=0.0125. Both retained a detected core at every sampled angle throughout their recorded histories. Maximum asymmetry was 0.106661266474 and 0.106661288252 respectively. Relative energy drift was 3.23e-11 and 1.01e-12 respectively. These are selected m=2 finite-time perturbations; twofold symmetry remains imposed by the initial data.

The original 16-angle and new 24-angle final states differ by at most 5.02e-11 at their eight coincident angles. This restricted comparison is not a continuum-convergence or unrestricted stability proof. The original interrupted 24-angle attempt remains excluded and is not replaced in its stop record.

## Physical background and binding

PHYSICAL_BACKGROUND.md reconstructs the physical potential from the implemented equations. The console ring sits on a charged rotating background. At lambda=25 an exact potential inequality gives E_relative >= sqrt(1.48)|Q_Phi| over the zero-charge vacuum. This excludes below-free-threshold binding in that stated charge sector, while allowing further questions about medium-supported or metastable structures.

The separately identified lambda=12 potential family admits a subthreshold energetic well. Compact spherical trial data and fixed-charge radial relaxation are exploratory calculations in that family, not changes to the console vortex equations.

## Spherical fixed-charge candidate

At Q_Phi=100000, Q_Sigma=0, lambda=12 and outer radius 60, radial relaxations reached residuals below 1e-7 on h=0.4 and h=0.2 grids. Energy/charge was 1.37736107558 and 1.37736653642, below the free Phi threshold sqrt(2). These are discrete spherical critical-point candidates. Relaxation iterations are not physical time. Domain convergence, additional grid refinement, full phase/angular stability, and comparison with measured particle properties remain unestablished. The calculation does not derive observed matter or establish a uniquely SEFI effect.

## Reproduction and integrity

Run test_extended.py and test_run_integrity.py with Python unittest. physical_potential_audit.py checks the force/potential reconstruction and exact algebraic identities. vacuum_trial.py reproduces compact trial energies. vacuum_spherical.py checks its variational gradient before reproducing both radial searches. These radial searches completed before the new pre-run marker mechanism was applied; no retrospective pre-run markers have been fabricated.

The angular run mechanism now writes a running status before computation, refuses reuse of a finished run identifier, and records interruption/stop status. snapshot.py admits only completed angular records that reached their requested duration. research/RESEARCH_DATA_POLICY.md records the standing exclusion policy. Existing historical paused-report PDFs retain their original scope; this note records the later results.
