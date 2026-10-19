# Verified background-response comparison

Jason Duran Dutton - 19 October 2026

This research checkpoint establishes an exact, reproducible difference between two specified scalar-field backgrounds: the gapless long-wavelength speed is 1 for the static global-Witten control and 1/sqrt(5) for the rotating console realization. It is a milestone in defining and checking the implemented model, not evidence of unique SEFI dynamics, nonlinear vortex stability, or empirical confirmation.

Read COMPARISON.md for the derivation, literature mapping and limitations. Run `python compare.py` with NumPy and Pillow to regenerate results.json and dispersion-comparison.png. Twelve mode cases check eigenvalues and physical-time integration of the linear Fourier system with time-step refinement. The console comparison is in SEFI-PY Evolution/sim/vortex_console at `/comparison`.

The original vortex-core survival question remains open. Next work must treat conserved charges consistently and test nonlinear, non-axisymmetric disturbances.
