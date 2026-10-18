# SEFI-PY repository audit and technical gap analysis

Jason D. Dutton | October 6, 2026 | Base commit 5b9ce48af8ed9f3a493bb6cf03c7227b194c02d6

## Audit scope and evidence

The live GitHub recursive tree contains 5,352 entries and is not truncated. All 256 repository-owned Python, JavaScript, TypeScript, Q#, and configuration sources in the filtered inventory were read, totaling 9,866 lines. Every inspected local source matches its pinned remote Git blob after allowing Windows CRLF normalization. Python files were statically parsed; this does not certify their runtime behavior. Dependencies, caches, rendered videos, binary archives, and the external gitlink are enumerated at tree level, not independently validated. `audit/source_inventory.csv` records every inspected path and definition; `audit/audit_record.json` records checks. `repository_inventory.json` retains owned paths and pinned hashes with dependency/cache/video paths excluded.

The repository has several overlapping trees, not one uniformly packaged Python library: `SEFI-PY/`, `SEFI-PY Evolution/`, a duplicated `SEFI QEC Stack/` subtree, JavaScript packages, and the separate `SEFI_QEC_Stack` gitlink. No new functionality is placed in a legacy import namespace.

## Existing components

| Area | Existing paths | Reuse assessment |
|---|---|---|
| Mathematical and identity modules | SEFI-PY/core/field_origin.py, field_authorship.py, field_sovereignty.py; Evolution/core/identity.py, sefi_space.py | Domain models, not validated quantum state spaces or channels. |
| Geometry | SEFI-PY/utils/geometry.py; sim/world.py; Evolution/geometry/curves.py, frenet.py, transforms.py | Vector operations and classical curves inspire reference-distance features. Avoid directly using legacy Frenet routines as quantum metrics. |
| Dynamical systems | SEFI-PY/sim/; Evolution/sim/defi_integrator.py, solver.py; core/time_dynamics.py | Classical evolution only. Preserve independently of quantum channels. |
| Observers and diagnostics | Evolution/research/calibration_workbench/workbench.py; SEFI-QEC/pauli_frame_check.py | Workbench already separates observations from truth and provides a useful workflow pattern. Pauli-frame checks compare simulated physical states to intended states; they are not experimentally measured observables. |
| Recurrence tools | SEFI-PY/core/time_dynamics.py | Exact rational recurrence, Floquet and RK4 diagnostics exist; no justified recurrence feature for this memoryless channel, so none is invented. |
| Phase tracking | time_dynamics.py and sefi-voxel/phase_compiler.py | Classical phase diagnostics and a compiler prototype; no quantum coherence guarantee. |
| Visualization | SEFI-PY/utils/plot.py, draw.py, draw3d.py; browser visualizers and Evolution/sim/vortex_console/ | Plotting patterns useful; scientific figures here use plain Matplotlib with independently defined metrics. |
| Testing | SEFI-PY/tests/, Evolution/tests/, calibration workbench tests, QEC Stack tests | Multiple isolated suites. Prior 33-test Time/core run is historical, not re-certified here. This addition has its own physical tests. |

## Missing components and misleading legacy names

| Required component | Present in some form? | What is still missing / minimal addition |
|---|---|---|
| Quantum state representation | Yes: SEFI-PY/core/quantum.py contains duplicate definitions of QuantumState, QuantumOperator and QuantumEvolution. | Later definitions overwrite earlier ones. States describe warp configurations; Euler stepping followed by normalization is not a verified unitary or general CPTP channel. Add a separate oscillator reference with exact overlaps and a finite-Fock matrix check. |
| Quantum channels | Bit-list noise and transport prototypes exist. | No validated bosonic random-displacement channel with unconditional fidelities. Implement a classical mixture of physical displacement unitaries and explicit CPTP erasure decoding. |
| Noise model | QEC Stack/quantum/noise.py includes bit flips; phase_flip is explicitly a placeholder. | Add a fully specified two-quadrature Gaussian displacement distribution and independent heterodyne shot noise. Do not treat bit-list phase metadata as dephasing. |
| Syndrome extraction | Evolution/qec/syndrome.py thresholds supplied displacement coordinates; legacy extractors exist. | No measurement instrument derives a bosonic stabilizer syndrome. The new benchmark uses auxiliary probe heterodyne outcomes, honestly labeled a displacement observation, not a stabilizer syndrome. |
| Baseline decoder | QEC Stack/quantum/decoder.py decodes 20 classical bits into Phi fields. | It is not a bosonic ML decoder. Add exact probe-likelihood ML and an independently implemented conventional Gaussian posterior. |
| Reliable quantum evaluation | Legacy quantum tests check amplitude changes/norms. | Add entangled-reference fidelity, six-state average, leakage, Kraus completeness, Choi positivity and cutoff convergence tests. |

Five legacy QEC Stack Python sources fail static parsing: azure_job.py, logical_register.py, registers.py, sovereignty_surface.py and syndrome_extractor.py. Paths and exact diagnostics are in the audit record. They are not imported, repaired, or silently reused in this benchmark. This is an isolated addition, not a certification of the legacy QEC stack.

## Minimal additions and implementation plan

1. Freeze the protocol and primary comparison before the simulation; enumerate all seeds and negative controls.
2. Add `bosonic.py`: oscillator states, displacement operators, exact code overlaps, unconditional metrics, CPTP erasure decoding.
3. Add `decoders.py`: measurement-only ML, geometry-distance regularization, independent conventional posterior, affine/nonlinear and randomized controls.
4. Add `run_benchmark.py` and `statistics.py`: paired held-out simulation, seed bootstrap, raw data, independent Gaussian quadrature check, timing and physical validation.
5. Add `test_benchmark.py`: independent checks of physics, inference and information boundaries.
6. Add `make_figures.py`, report sources, reproducibility notes, raw/aggregate results and publication-quality PNG/SVG figures.

No modification of the legacy core is necessary. The new directory is `SEFI-PY/benchmarks/displacement_recovery/`; it can run standalone from a repository checkout. Reuse is conceptual and architectural: observation/truth separation and reference-geometry distances. No SEFI equation is asserted to generate a quantum channel.

## License and practical limits

The existing root license is SEFI-PY Research License v1.0, allowing noncommercial research with attribution; it is not a permissive commercial open-source license. This task does not change that license. No cloud platform, hardware access, or paid tool is needed to reproduce the benchmark. A new grant-funded licensing plan remains a separate organizational decision.
