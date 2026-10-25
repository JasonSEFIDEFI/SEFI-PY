# Execution and preservation record

All runs are exploratory, dimensionless and adaptive, not preregistered.

- frequency_scan_initial.py generated scan_X32. Its larger frequency steps sometimes exhausted nodes or converged to the trivial vacuum, especially near lambda=16. Rejected attempts remain in its JSON and NPZ outputs.
- frequency_scan.py changes seed amplitude scaling: it is applied only for beta<0.08. Applying it indiscriminately to the strong-amplitude branch was an avoidable numerical seed problem.
- recover_critical.py used smaller frequency steps at lambda=15.5,15.75,15.9 and recovered nontrivial branches. This does not alter the governing equations.
- verify_turn_initial.py obtained the lambda=12 minimum and pair and some binding crossings, then stopped when a lambda=15 solve returned the trivial vacuum. Its partial records are retained under verification_initial.
- verify_turn.py uses small-step continuation and a tighter binding bracket. Its completed records are under verification_refined.
- The first paired-Hessian pass completed the high-frequency checks, then stopped after an input path was relocated. The complete rerun for both profiles is paired_hessian_refined.json; the earlier high-frequency-only record is retained.
- critical_limit_initial.py reached a failed large step to beta=0.05. Its accepted preceding points are retained in critical_limit_results.json; the failed solve had status 1, maximum residual about 0.33, and was not accepted as a reduced-model solution.
- critical_limit.py uses smaller steps through that region and completes the fold and binding calculations in critical_limit_refined.json.

Original public studies and prior checkpoints are not overwritten. The verification/ high-frequency profile is a duplicate of the copy retained in verification_initial and is omitted from the publication package to avoid duplicating an identical input.

