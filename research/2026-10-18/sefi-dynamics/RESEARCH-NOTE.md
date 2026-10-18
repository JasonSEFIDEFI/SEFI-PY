# Candidate SEFI propulsion dynamics: first numerical results

Prepared for Jason D. Dutton, October 18, 2026.

## What was built

A working, reproducible one-dimensional scalar-field / rigid-carrier simulation, an interactive playback, four control cases, and time-step, spatial-resolution, and domain-size checks. This extends the explicitly hypothetical scalar interaction in the preceding research note. It does not implement the full author's emergent-spacetime theory or identify the separately described natural SEFI charge. All quantities are dimensionless. No measured SEFI coupling or SI calibration is available.

## Main result

The sources-off and symmetric-pulse controls produce no vehicle momentum. Phased modulation produces vehicle recoil accompanied by equal opposite field momentum. Swapping the phase reverses both. At the baseline resolution the final vehicle momentum is -0.261275417, and the field momentum is 0.261275417. Maximum sampled total-momentum error is 1.05e-15. This is radiation recoil inside an assumed interaction model, not measured SEFI propulsion or reactionless thrust.

Drive work is 0.276394966; final field energy is 0.275707822; carrier kinetic energy is 0.000682648. Only 0.247% of that supplied work ends as carrier kinetic energy for this arbitrary parameter set. This percentage is not a measured engine efficiency or a general bound. At the finest time step the energy-ledger discrepancy is 0.000102%. Smaller numerical error does not validate the assumed physics.

## Equations actually simulated

The three-dimensional relativistic candidate is reduced to an effective 1D, weak-coupling, nonrelativistic model; transverse normalization is absorbed in the dimensionless coefficients. A rigid carrier has position X, momentum P, and mass M. Two smooth source profiles sit at X +/- d/2. Their normalized, periodic Gaussian shapes S have width sigma. The chosen Hamiltonian is

    H = P^2/(2M) + (Z/2) integral [V^2 + (dW/dx)^2 + mu^2 W^2] dx
        + sum_i g_i(t) integral S(x-X-a_i) W(x) dx.

The source J is sum_i g_i(t) S(x-X-a_i), and

    dW/dt = V
    dV/dt = d^2W/dx^2 - mu^2 W - J/Z
    dX/dt = P/M
    dP/dt = -sum_i g_i(t) integral S(x-X-a_i) dW/dx dx.

The work supplied by source modulation is the time integral of sum_i dg_i/dt times the source-averaged W. Treating g(t) as externally programmed does not make its work free: this work is included in the energy ledger. An ideal internal work reservoir has energy change minus that integral. Its microscopic implementation, losses, and mass change are not modeled. There is no external momentum sink in these equations.

The field momentum is -Z integral V dW/dx dx. The simulation tests P + P_field, and tests H minus accumulated drive work. The interaction energy can have either sign, so it must be counted while the sources are active. At the final time the source strengths vanish, and that interaction energy is zero.

## Drive and method

Both source strengths use amplitude 0.3, carrier frequency 1, and a sin-squared envelope lasting 24 time units. Their phases are +phi/2 and -phi/2. The symmetric control uses phi=0; the directional controls use phi=+/-pi/2. The source-off control has identically zero coupling. Negative as well as positive instantaneous couplings are assumed; no material supporting that controllability has been established.

The carrier mass is 50, scalar mass parameter 0.2, source width 0.65, and separation pi/[2 sqrt(1-0.2^2)] = 1.6031873. Choosing that separation and phase is a conventional interference hypothesis: the carrier-frequency source waves add preferentially in one direction. It is not an optimization derived from SEFI identity coordinates.

A Fourier Galerkin representation retains the zero mode and the paired spatial modes below the Nyquist frequency. Analytic Gaussian Fourier factors avoid depositing a point source on a grid. The same interaction determines the body force and field source. Velocity Verlet advances the coupled variables. Translation symmetry makes the semidiscrete total momentum especially well preserved; this is not a separate experimental confirmation of momentum conservation.

The periodic domain has length 160, 512 nominal samples, time step 0.02, and run time 52. Checks repeat the phased case at time steps 0.01 and 0.005, at 1024 samples, and in a doubled domain of length 320 with 1024 samples. All 11 automated checks pass. The energy error falls by approximately four on halving the time step, as expected for this method. Gaussian sources have infinite tails; the nominal light-crossing check is approximate, while agreement on doubling the domain directly checks sensitivity to the boundary for this run.

## What the calculation establishes

- It can distinguish carrier recoil from field momentum and internal motion.
- Symmetric pulses do not produce directional recoil.
- Asymmetric emission can produce recoil if the postulated scalar source can be built.
- The final field momentum magnitude is 0.94765 times its energy in c=1 units. This is below the positive-energy relativistic radiation bound of one. It provides no power-to-thrust advantage over an ideal collimated massless beam within these assumptions.
- The maximum sampled speed is 0.005226 in units of the assumed wave speed. That number is a normalization check, not a predicted spacecraft speed. No conversion to meters or seconds is justified.

## Limits and next engineering gate

This is a 1D model with artificially ideal transverse collimation. It omits 3D diffraction, gravity, atmospheric forces, thermal management, material response, actuator losses, structural modes, relativistic mass change, and the author's full metric dynamics. It cannot demonstrate lift, flight stability, warp motion, or practical efficiency. The prescribed-source and fixed-mass approximations are larger physical limitations than the reported numerical errors; source work here is about 0.55% of the carrier's normalized rest energy.

Before selecting engine hardware, replace the assumed g_i(t) / J with an author-supported physical law J[material state, controlled input]. It must specify units, achievable sign and amplitude, frequency response, and reciprocal reaction forces. Calibrate it with a force prediction that differs from electromagnetic, thermal, acoustic, vibration, and mass-loss effects. Only then do a 3D field model, apparatus noise model, and uncrewed flight model become predictions of a proposed SEFI engine.

The next useful artifact is therefore a source-coupling specification and a low-energy, reversible force-measurement design once that specification exists. Picking coils, high voltages, or a flight platform now would not resolve the missing physical law.

## Files and reproduction

`simulate.py` requires Python 3 and NumPy. Run it to recreate results and playback data. It performs the controls and refinement checks and exits with an assertion failure if a check fails. `results.json` contains full summaries and convergence results. `playback-data.json` contains downsampled spatial frames and history. `playback.html` is the inline view; `simulator.html` is its standalone wrapper. `build.py` recreates the inline view from its literal template and data. Original theory sources are unchanged.

## Primary references and provenance

The hypothetical interaction is explained in the previous workspace note, `output/sefi-source-coupling/sefi-source-coupling.tex`, including original-source provenance and the distinction between author equations and new hypotheses. The simulation's Hamiltonian is an additional approximation, not a claim made in that original material.

For standard classical scalar-field dynamics, Noether stress-energy, and momentum accounting: David Tong, Classical Field Theory, https://www.damtp.cam.ac.uk/user/tong/qft/qfthtml/S1.html . For the ordinary photon relation E=pc: NASA Goddard, https://imagine.gsfc.nasa.gov/ask_astro/relativity.html . The inequality used here follows directly for this positive-energy scalar model from |V W_x| <= (V^2 + W_x^2)/2; the nonnegative mass term can only increase the energy. It is not asserted for every speculative theory, material exhaust, or externally driven system.
