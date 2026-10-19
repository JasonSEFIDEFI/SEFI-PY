# Experimental nonlinear angular extension - paused checkpoint

Jason Duran Dutton | 19 October 2026

The user stopped this batch early. It did not exhaust the authorized 30-minute allowance. See STOP_RECORD.json for complete, interrupted and unrun cases. No simulation is scheduled to resume automatically.

The separate cylindrical prototype extends the preserved lambda=25, N=0 equations into an angularly periodic, nonlinear r/theta/Z evolution. It uses finite differences in r,Z and a Fourier angular Laplacian, with axis regularity enforced by a single angularly averaged axis value. It is not an Authorship-dependent model. The main console's physics.py, ring_operator.py and rest_search.py are unchanged.

Four prototype tests pass: reduction to the axisymmetric equations; a harmonic r^2 exp(2i theta) angular-Laplacian check; directional agreement with the prior angular tangent equations; and differential discrete energy/carrier-charge checks. These are consistency checks, not empirical validation.

Completed pilot control: 8 angles, T=8, no imposed disturbance; all sampled meridional angles retain a detected core. Pilot m=2 disturbance: 8 angles, T=8. Main m=2 disturbance: 16 angles, T=12; all sampled angles retain a detected core, maximum angular asymmetry gain 1.2084, relative energy drift 3.23e-11 and carrier-charge drift 2.77e-11.

The 24-angle run was interrupted without saving its in-progress state. Angular and temporal convergence have therefore NOT been established. Pilot and main runs have different final times and cannot substitute for a matched refinement comparison. The selected m=2 initial condition retains a discrete symmetry and does not test arbitrary 3D disturbances. Nonlinear collocation products can alias; more angles, spectral diagnostics and matched runs are needed. Fixed boundaries, finite time, incomplete spatial convergence and absence of Authorship coupling remain substantive limitations.

Report.pdf records both this partial batch and the earlier completed survival study. The console receives a read-only evidence snapshot; experimental solvers do not become selectable production models.

Reproduction: install NumPy; run `python -m unittest test_extended -v`. Invoke `cylindrical3d.run_case` with explicit parameters for an individual case. `batch.py` preserves the dated plan and original cutoff; do not assume it is a currently scheduled computation. `snapshot.py` regenerates the read-only evidence snapshot from saved completed cases. `record_stop.py` verifies protected equations and records this stop checkpoint. Pillow and ReportLab are needed only for figures/reporting.
