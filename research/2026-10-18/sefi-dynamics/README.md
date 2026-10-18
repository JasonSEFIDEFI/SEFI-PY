# Candidate field-recoil simulator

This package contains a working 1D scalar-field simulation and an interactive **playback of computed results**. The controls choose between already-computed scenarios and time slices; they do not rerun the field solver. The physics is an explicitly hypothetical interaction, not a validated SEFI engine or a flight simulation.

Open `simulator.html` in a browser. It uses an approved public CDN for the D3 chart library, so the interactive charts require an internet connection. Select a run and move the time slider to inspect its field and momentum history. The shown quantities are dimensionless, not measured vehicle performance.

## Reproduce the numerical results

Install Python 3 and NumPy (`python -m pip install numpy`), then run:

```
python simulate.py
```

The program runs four scenarios and numerical refinement checks, writes `results.json` and `playback-data.json`, and raises an error if any of its 11 checks fails. Runtime on the author's workspace was approximately 30 seconds. The simulation does not need access to the V drive.

`RESEARCH-NOTE.md` gives the exact assumptions, equations, parameter values, results, and limitations. `MODEL-ORIGIN.tex` is the preceding mathematical note explaining the newly proposed interaction and the distinction from the author's source theory.

`python build.py` rebuilds the inline playback fragment and research note. The included standalone browser wrapper was generated with the workspace's visualization tooling. The automated browser-check results are recorded in `ui-checks.json`; the browser script is workspace-specific.

No source material was changed, and nothing was submitted or published. A source-to-material coupling and experimental calibration are still required before this can predict physical thrust.
