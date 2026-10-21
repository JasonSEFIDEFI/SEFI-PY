# Photonic theory-to-experiment bridge

Local sources inspected October 21, 2026:

- `V:/V Drive and GitHub Content/Journal Submissions/Photonics Submission/manuscript.tex`, particularly stability geometry, correction mapping and photonic implementation.
- Header of `V:/V Drive and GitHub Content/Journal Submissions/FULL SUBMISSION/manuscript.tex`; these local files have differing lengths. The exact revision under review has not been verified.
- Earlier consolidated audit `output/unified-jmp/unified-manuscript.tex`, particularly its explicit limitations on photonic interpretation and projection.
- SEFI-PY's time-dynamics reproduction guide and the displacement-recovery benchmark.

User-reported submission history: photonics manuscript transferred from AVS Quantum Science to AIP Advances and under review there. No APL Photonics submission is asserted.

## What follows and what does not

The photonics source proposes an optical-mode configuration space, stability set S(x)=0, geometric departures, and correction back toward a reference set. It does not uniquely specify a measured state x, explicit constraint S, metric/domain for projection, quantum instrument implementing the syndrome, or recovery channel with a stated noise model. Calling curvature/torsion a syndrome does not supply that instrument. A mode-space trajectory need not be a spacetime photon worldline; definitions and units must be selected explicitly, not silently transferred from timelike proper-time formulas.

The intended retention/restoration principle can motivate an engineering objective. The intervention gate gives that objective a physical overlap score in a separately defined oscillator model. Its extra ingredient is a calibrated actuator-disturbance model. That ingredient is not derived from the original curvature/torsion geometry and is explicitly an extension.

## A genuinely SEFI-specific prediction would require

1. A specified optical state and measurable operational meaning for the proposed identity variables, including units and calibration.
2. An evolution or constraint equation giving a numerical optical prediction, not simply a label for ideal operation.
3. A valid physical measurement/recovery map and proof it respects positivity, normalization and information access. If conditional filtering is used, success probability and failed trials must be included.
4. A quantitative difference from the corresponding conventional optical/control model under the same independently measured parameters. Free coefficients cannot be adjusted after observing the target signal.
5. A small reproducible test separating that difference from birefringence, ordinary geometric phase, attenuation, detector errors and familiar feedback effects.

No rule is introduced merely to ensure a large favorable benchmark. A usable engineering improvement can be valuable without being a new physical effect, but its provenance must remain clear.

## Established geometric-optics overlap

Curvature/torsion effects on polarization and geometric phase in helical fibers already have quantitative treatments. For example, Band, Kuzmenko and Avishai, *Geometric phases in optics: Polarization of light propagating in helical optical fibers*, Phys. Rev. A 111, 033530 (2025), https://journals.aps.org/pra/abstract/10.1103/PhysRevA.111.033530 ; accessible author preprint https://arxiv.org/abs/2407.19425 . The model includes stress-induced birefringence and warns that geometric phase is not automatically topologically protected. Reproducing such a phase would not by itself be a SEFI-specific discovery.

No submitted manuscript was edited; no experimental invitation, journal submission or equipment purchase was made.
