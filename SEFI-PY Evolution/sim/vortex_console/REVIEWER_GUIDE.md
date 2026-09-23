# Reviewing the vortex research

This is an exploratory classical field model. Its motivating primordial-unity hypothesis is not an established premise of physics. A numerical vortex is not yet an identified physical particle. Review should distinguish hypotheses, derivations, finite-grid observations, and unresolved claims.

## Using GitHub Copilot after publication

Open the published repository in GitHub, open Copilot Chat, and include the repository or specific files as context. Reviewers use their own available Copilot access; publishing files does not share the author's personal Copilot account, chat history, subscription, or private manuscripts. The instructions file guides supported Copilot experiences; support varies by interface and instructions cannot guarantee an accurate answer. Check response references to see which files were used.

The questions below are also accessible from the console's **Review** menu. Anyone can read this guide without Copilot and use the same questions for a human review. No hosted chatbot or external AI connection is embedded in the console.

Official documentation: [asking questions on GitHub](https://docs.github.com/en/copilot/how-tos/copilot-on-github/chat-with-copilot/chat-in-github), [repository instructions](https://docs.github.com/en/copilot/how-tos/copilot-on-github/customize-copilot/add-custom-instructions/add-repository-instructions), [instruction support](https://docs.github.com/en/copilot/reference/custom-instructions-support).

## Common instruction to append

Use this repository as evidence. Cite file paths and relevant equations or tests. Separate established physics, model assumptions, derived results, numerical observations, and speculation. Do not treat SEFI, DEFI, or GWFM as established physics premises. Distinguish stored results from calculations you actually ran. If a source or test is missing, say so. Identify uncertainties and constructive next tests rather than assuming the intended conclusion.

## Questions by reviewer perspective

### General reader or collaborator

1. What does this model actually establish today? Explain it without specialist vocabulary, then give the mathematical version. What does the rendered surface represent?
2. Which controls change the observation, which change the equations, and which change only an initial guess? Demonstrate with a reproducible example.
3. How far are these results from deriving physical matter? Identify the missing steps without promising that the research will succeed.

### Theoretical physicist

1. Derive both field equations from the stated Lagrangian, with signs, normalizations, units, background, and boundary conditions explicit. Does the traveling ansatz reproduce the implemented equations?
2. Is the carrier charge a model Noether charge or a derived electromagnetic charge? What further structure would make an identification with observed particles justified?
3. What would be needed to establish mass, spin, quantum statistics, interactions, and experimentally measured scales? Which properties are absent from the current model?
4. Does the effective phase metric apply only to one mode? What calculation would test whether all sectors share the same propagation geometry? Is any gravitational field equation derived?

### Mathematician

1. Audit the rest-state obstruction, including finite energy, background, time dependence, boundary terms, and scaling assumptions. Identify counterexamples outside its scope.
2. Does the proposed alternative parameter window remove an obstruction or prove existence? Formulate the next meaningful existence or nonexistence question precisely.
3. Distinguish a small discrete residual from a continuum solution and from nonlinear orbital stability. What compactness, error estimates, or spectral analysis are missing?
4. Explain the difference between phase winding around the vortex core and carrier winding around the symmetry axis. Which topological properties survive the chosen boundaries?

### Experimental physicist

1. Which computed quantities could map to measured observables, and what independent calibration of dimensionless scales is required?
2. Is there a presently derived falsifiable prediction that distinguishes this model from established physics? If not, what is the smallest calculation needed to obtain one?
3. What observations could reject the model? Separate adjustable parameter fits from predictions made before fitting.

### Numerical analyst

1. Reproduce the baseline residual, charge, and ring radius. Compare grid refinement and domain enlargement separately. Quantify the remaining boundary sensitivity.
2. Repeat an unperturbed control and time-step refinement with a nonzero perturbation. Report actual errors and conservation drift, not just successful completion.
3. Compare raw carrier distance with global-phase-aligned distance. Could apparent deformation instead be phase drift?
4. Which unstable modes cannot be represented by axisymmetric evolution? What would a full three-dimensional perturbation study add?
5. Verify checkpoint continuation against uninterrupted evolution using the complete complex fields, velocities, time, and original charge normalization.

### Software and reproducibility reviewer

1. Trace parameter input through validation, nonlinear solution, integration, diagnostics, stored frames, and rendering. Are any visual changes substituted for physical calculations?
2. Check failed-solve labeling, cancellation, pause/resume, parameter sweeps, storage limits, and restoration after server restart. Provide exact reproduction steps for failures.
3. Check that the console runs independently of the SEFI-PY engine, uses only loopback networking, and does not expose arbitrary files or perform external AI calls.
4. Can a reviewer reproduce the claimed test results from the documented environment and supplied baseline data? Distinguish executed tests from proposed future tests.

## Evidence map

| Question | Starting evidence |
|---|---|
| Equations and numerical controls | `physics.py`, `ring_operator.py`, console Math menu |
| Saved profile provenance | `data/charged_ring_*.json` and matching `.npz` |
| Run lifecycle and checkpoint data | `server.py`, local `runs/<id>/settings.json`, `result.json`, `checkpoint.npz` |
| Visualization versus physics | `static/app.js`, `static/index.html` |
| Executed local verification | `TEST_REPORT.md`, `test_console.py` |
| Larger research claims | The repository's existing `research/physical_matter_2026_09_23` checkpoint and manuscript; cite the exact file actually consulted |

Run files are local and not automatically published. Ask for the relevant run artifact when it is absent; do not infer results from an image alone.
