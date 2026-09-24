# Local verification — 23 September 2026

## Console 0.2 checks — 24 September 2026

The six original numerical/API tests passed during the combined regression run. The new integer-grid test exposed integer truncation of the axis quadrature weight when Search was constructed directly with integer h=1. The constructor now makes grid coordinates floating point; all four fixed-charge tests subsequently passed in 7.216 seconds. The previously reported research runs used fractional h and were unaffected. The new tests check energy-gradient consistency, actual energy descent, input bounds, exact checkpoint continuation, model labeling, and the saved refined candidate. Existing physical-evolution equations are unchanged.

The portable two-charge audit ran successfully and verified charge compensation and its energy identity for five support radii. The live API returned the expected refined candidate, and a short search was submitted through HTTP. JavaScript syntax checks passed. The updated browser could not be visually verified because browser automation timed out twice; visual QA for the Charge tab remains outstanding. The earlier browser checks below apply to console 0.1, not proof of the new tab's visual quality.

Environment: Windows; Python 3.12.14; NumPy 2.3.5; BLAS thread count 2 for tests. Browser checks used the Codex in-app browser. Source version: coupled-ring-axisymmetric-v0.1.

## Numerical and API checks executed

The initial four-test suite passed in 40.436 seconds. It covered parameter validation, nonlinear parameter response, short evolution/control/refinement, and API/run/checkpoint behavior. Additional failure-path and model-version checks are in the same reproducible suite.

| Check | Observed result |
|---|---|
| Saved refined-profile stationary RHS | Maximum absolute residual 3.8628211740e-11 |
| Unperturbed evolution through T=1 | Original-field weighted distance 9.8453443e-11; relative charge drift 1.2450818e-16 |
| Perturbed evolution through T=1, ε=.001 | Core winding 1; relative charge drift 3.6617038e-13 |
| Time-step refinement .025 versus .0125 through T=1 | Maximum difference across complex fields and velocities 9.6126419e-10 |
| Changed travel coefficient c=.58 | Newton residual 5.7032157e-13; ring radius 5.6652955, compared with saved baseline 5.3574045 |
| Checkpoint at T=.5, continued to T=1 | Maximum state difference from uninterrupted run 5.5512448e-17; original charge normalization retained |
| Independent ε sweep: 0 and .002 | Both completed; unperturbed control remained near stationary and perturbed field changed |
| Lifecycle | Pause, resume, stop, saved checkpoint verified through local HTTP API |
| Rejected requests | Unsupported evolution λ/N, malformed sweep, missing local-request header, invalid checkpoint path |

Failed-solver behavior is tested with an injected nonconvergence result: it must preserve a labeled trial frame and never create an evolution. This is a software guard test, not a new physical rest-state search. A synthetic incompatible checkpoint tests model-version rejection.

## Browser checks

The saved profile loaded with the expected charge, winding, and residual. The top camera control changed elevation from 28° to 89°, and visual inspection confirmed the top view. Reset returned the default view. A T=1 evolution executed through the UI, completed with five saved frames, displayed charge drift 3.662e-13, and exposed checkpoint/settings download links. Replay entered playback mode. No warning or error messages were recorded in the browser log during that run.

The Review menu exposes six reviewer perspectives. The interface does not claim to be a connected Copilot chat. Additional human usability review is welcome, particularly transparency choices and small-screen layout. Not every browser, camera gesture, export path, parameter combination, storage exhaustion case, or operating system has been tested.

## Scope

These checks establish the stated local implementation behavior and short-time numerical consistency. They do not establish continuum convergence, long-time nonlinear stability, full three-dimensional stability, observed physical matter, universal matter coupling, or gravitational field equations. The earlier research checkpoint supplies longer-time and domain studies; those are separate evidence and were not all rerun for this console release.

Reproduce with `python test_console.py`. The tests use a temporary run directory and do not modify the baseline data or existing engine.
