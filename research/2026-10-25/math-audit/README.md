# Reproduce the mathematical audit

Install the versions in requirements.lock.txt in a separate environment.
The prior critical-binding package must remain alongside this output package, or in the repository at research/2026-10-23/critical-binding.

Run from this folder:

    python continue_math.py
    python refine_landmarks.py
    python second_alpha.py
    python stability.py
    python foundations_floquet.py
    python -m pytest test_audit.py -q

The committed selected test logs for QEC, geometry, nuclear regression and Voxel refer to their separate source packages, not this folder alone. Re-run those packages' test paths to regenerate the corresponding logs.

Use OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1 to reproduce the captured threading configuration. Timings are machine observations, not performance guarantees.
The scripts overwrite their generated outputs. Copy the package before rerunning if you want to keep the captured reference results.
summarize.py requires all result files and the five selected-suite logs and generates the report, figures and integrity manifest. Failure logs describe earlier attempts and are intentionally retained.

The mathematical claims and their assumptions are in MATH_CONTINUED.md. The scope and observed outcomes are in REPORT.md.
