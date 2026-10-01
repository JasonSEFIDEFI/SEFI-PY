# Sparse moving-profile solver: measured result

Tested September 30, 2026 on Windows-11-10.0.26200-SP0, Python 3.12.14, NumPy 2.3.5, SciPy 1.18.1, with OPENBLAS_NUM_THREADS=2 and OMP_NUM_THREADS=2.

## Matched computation

Saved charged-ring reference, n=32, L=16, c=0.58, nu=0.3, lambda=25, N=0, seed=5.36. Each backend starts independently from the same reference. Three runs per backend alternate order; elapsed times include matrix construction, the Newton solve, and the first sparse dependency import. Median timings describe repeated solves, not cold startup. No regression or competing benchmark process ran during this final measurement. Raw samples and full settings are in SPARSE_BENCHMARK.json.

| Measurement | Dense | Sparse |
|---|---:|---:|
| Median total profile solve (seconds) | 16.904947 | 0.245975 |
| Stored linear operator (bytes) | 73932800 | 236828 |
| Final maximum equation residual | 5.70322e-13 | 5.69655e-13 |
| Newton updates | 5 | 5 |

Measured median speedup: **68.73x**. Largest absolute difference across all u, w, s grid values: **4.996e-16**. Linear-operator storage fell by **312.18x**. Operator storage excludes Jacobian factorization, arrays, allocator overhead and process memory; it is not a peak-memory measurement. Speed depends on hardware, library builds, grid, convergence and system load; this one converged profile does not establish a universal speedup.

## Implementation and checks

The original finite-difference stencil is assembled directly as sparse triplets. Newton adds the same diagonal and interfield derivatives and solves the resulting sparse system. Winding constraints replace the same rows. The original dense implementation remains selectable and is the automatic fallback without SciPy. Equations, boundaries, convergence tolerances, backtracking and the checkpoint model version are retained.

All 18 checks passed: 15 existing numerical/API/checkpoint/fixed-charge/warp tests and 3 added equivalence/fallback/cancellation tests. Stencil entries and boundary source vectors match exactly at multiple grids and travel coefficients. Solved fields agree within rtol=1e-7, atol=1e-8 for n=16 and 32, including N=1 boundary enforcement; these test cases compare solver outcomes rather than claiming every trial converges. The n=32, N=0 benchmark converges with both engines.

Reproduce from this directory:

```text
python -m pip install -r requirements-fast.txt
python -m unittest test_console test_rest_console test_warp_vortex test_sparse_solver -v
python benchmark_profile_solver.py
```

SciPy API references: [sparse direct solve](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.spsolve.html), [triplet assembly](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html).

## Meaning for the theory

This is a computational improvement in the existing classical coupled-vortex model. Faster profile searches reduce the cost of parameter exploration and comparison. They do not establish physical matter, continuum existence, nonlinear stability, a propulsion effect, or new gravitational dynamics. Fixed-charge relaxation and physical time evolution were regression tested but are not accelerated by this change.
