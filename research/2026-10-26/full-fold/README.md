# Reproduce this session

Captured Python version: 3.14.7. Use requirements.lock.txt and retain the prior audit's dependency folders:

- output layout: this package beside sefi-math-audit-2026-10-25 and sefi-critical-binding-2026-10-23.
- repository layout: research/2026-10-26/full-fold, research/2026-10-25/math-audit, research/2026-10-23/critical-binding.

Set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1.

Run in sequence:

    python session.py
    python radial.py
    python sensitivity.py
    python -m pytest test_math.py -q
    python summarize.py

Capture the test output as tests.log before running summarize.py. Scripts overwrite their output files; use a copy to retain the saved reference results.
No hardware, email, external publishing or device action occurs inside these experiment scripts.
All quantities are dimensionless. See MATH.md for derivations and REPORT.md for observed results and limits.
