# Time update validation — October 6, 2026

Author: Jason D. Dutton.

## Results

- 12 new Time dynamics tests passed.
- Full SEFI-PY suite: 33 tests passed (Python 3.12, pytest 9.1.1; the existing engine tests also used NumPy and Matplotlib).
- The standalone demonstration reproduces the resonant monodromy matrix to displayed precision, Floquet growth rate 0.035314022950264, stable comparison multipliers, nonlinear barrier 0.5, and exact rational recurrence certificate at 10*pi.
- The new dynamics module and demonstration require only the standard library. No new runtime dependency was added to the engine.

Tests cover the analytic spectrum and orthonormal eigenbasis across random signed couplings, agreement between exact linear trajectories and independently integrated equations below/at/above the stability threshold, completed-square potential identity, stationary saddle and center, the prescribed-work energy identity, RK4 convergence, determinant conservation, harmonic phase lifting across multiple cycles, sampling ambiguity rejection, exact rational recurrence certificates, nonperiodic observable examples, and DNA energy derivatives.

Initial broad testing exposed an existing missing import in `tests/test_core.py`: it instantiated `WarpExpression` while importing `WarpEngine`. The import was corrected, without changing engine behavior. All tests subsequently passed.

Numerical checks establish implementation consistency with the stated mathematical model. They do not establish experimental validity, global coercivity of arbitrary coupling functions, physical clock calibration, or attracting synchronization. The manuscript is described as submitted, not accepted.
