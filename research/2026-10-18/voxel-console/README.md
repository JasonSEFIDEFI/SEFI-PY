# SEFI Voxel Machine Console

An offline, dependency-free interactive local research console. **SIMULATION ONLY. No connected hardware, measured positions, calibrated pressure or physical levitation evidence.**

## Run

Requires Node.js 20 or later. In this folder run `npm start`, then open http://127.0.0.1:8765. Run `npm test` for the automated model and hardware-boundary tests. No installation step or network dependency is needed. `PORT` can select another port. The server binds only to loopback. Stop the server with Ctrl+C.

## Explore

1. Press Run simulation. Drag the workspace to orbit, scroll to zoom, double-click to reset the view.
2. Select circle, square, line or helix, edit radius/period, and Apply configuration. Then Run. XYZ offsets move the trajectory center. Helix traverses two turns upwards and returns along the same path for continuous looping.
3. Change array geometry, drive settings and particle parameters. Apply validates simulation boundaries and resets telemetry. Field slice shows relative XY phasor magnitude at commanded Z; no pressure units. Export channel phases for inspection, not as a validated hardware recipe.
4. Record a VM-EXP ID and export JSON/CSV. Recording occurs at 10 Hz simulated time; physics advances at 100 Hz. Pausing stops simulated time. Notebook and inventory persist in this browser origin. Export before clearing data. The latest 30 completed experiments persist; active recordings remain in memory until finished. Each recording is capped at 36,000 samples.
5. Run VM-EXP-001 to evaluate a static hold after configurable settling. PASS/FAIL applies only to synthetic data. The physical protocol remains NOT TESTED.
6. SIM E-STOP latches simulation off. Reset re-arms it; it cannot cut physical power.

## Models and evidence

Array coordinates are millimeters, frequency Hz, phase degrees, amplitude normalized 0–1, particle diameter mm, density kg/m³. Array planes lie at ±gap/2 along Z. The point-source field uses c=343 m/s and distance-based focusing phases with a configurable opposed-plane offset. It computes a normalized phasor magnitude, without directivity, scattering, pressure calibration or near-field validity. A focus phasor is not a verified stable acoustic trap.

Synthetic particle motion is an independent arbitrary damped PD response: a=Kp·R·(command−position)−(Kd+damping)·velocity. R depends on amplitude, diameter, density and array layout. It is **not** derived from radiation forces. Frequency and phase affect only the field. Gravity, Gor’kov potential, real trapping, thermal effects, measurement noise, camera calibration and physical stability are absent. Numeric workspace clamping prevents simulator divergence; it cannot certify physical containment. Instant/RMS/max error compare synthetic position against the speed-limited command, not the raw requested trajectory. RMS and maximum cover all physics steps since reset.

DERIVED / SIMULATED / COMPUTED / BENCH-TESTED / MEASURED / REPRODUCED / CONJECTURED are distinct evidence labels. This app generates only SIMULATED and COMPUTED data. Physical inventory status is a user assertion, initially PLANNED throughout. The SEFI research section stores a hypothesis; it never changes the baseline controller and makes no performance claim.

## Integration boundary

`model.js` exports SimulationAdapter and an intentionally unimplemented HardwareAdapter. Hardware connect/command reject; measure returns null and disconnected; stop returns no physical acknowledgment. Before enabling outputs, implement a separate adapter with capability discovery, calibrated units, channel mapping, acknowledgment/timeouts, independent physical cutoff, driver protection and timestamped camera measurements. Keep SIMULATED provenance separate from MEASURED evidence. No serial/camera permissions are requested by this app.

Progression: simulation → bench-validated first levitation → hardware control → calibrated camera feedback → measured closed-loop positioning. Neither the digital twin nor a simulated PASS advances physical readiness.

## Files

`index.html`, `style.css`, `app.js`: user interface, projection/orbit rendering, notebook and exports. `model.js`: geometry, phase calculations, trajectories, synthetic dynamics, criteria and adapters. `server.js`: local static server. `model.test.js`: tests for bounds, trajectories, convergence, metrics, criteria and disconnected hardware. No synced project source files are modified.
