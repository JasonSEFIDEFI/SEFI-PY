# SEFI Double Slit Explorer

A five-chapter educational application: begin with the double-slit experiment, then examine different ways of describing it. Built with React, TypeScript, Vite, Tailwind CSS, Framer Motion, and Canvas 2D.

**This is a separate educational app. It does not modify the SEFI-PY engine, SEFI-QEC, or the Vortex Research Console.**

## Run locally

Install Node.js 24 LTS, open a terminal in this directory, then run:

```sh
npm install
npm run dev
```

Open http://127.0.0.1:5173. To stop the server, press Ctrl+C in its terminal. See [INSTALL.md](INSTALL.md) for production and hosting instructions.

```sh
npm test
npm run build
npm run preview
```

The production output is `dist/`. No server, API key, account, or external AI service is required. The app itself makes no analytics calls. External reference links open only when selected.

## The experience

1. **Quantum mechanics:** interference amplitudes, diffraction envelope, and individual detection events. Slit width and separation change a calculated probability distribution.
2. **Quantum field theory:** an educational field analogy, with interactions and a path marker. It is not a numerical QFT solver.
3. **GWFM:** an explicitly proposed path-bundle interpretation. Path count, spread, and overlap are drawing controls, not measured physical quantities.
4. **SEFI:** Origin, Authorship, Sovereignty, and Warp as interactive identity concepts. These are proposed interpretive labels, not established quantum observables.
5. **DEFI:** an illustrative geometry–identity coupling, with visible input and response. It is not a derived physical coupling law.

The shared observation slider updates all five explanation cards. The Advanced Research View is closed initially, keeping equations out of the first-time experience. A final research section frames geometry → particles → matter as open questions.

Controls include pause/play, a batch of 1,000 detections, clear screen, image export, slit width, slit separation, marker enablement, and marker interaction. Physical changes restart the deterministic sample; switching chapters retains it. Animation initially pauses when the system requests reduced motion. Native controls support keyboard use and visible focus.

## Scientific boundary

Every chapter uses the **same standard quantum reference distribution** for screen detections. The geometric animations do not independently derive that distribution. The information slider is an idealized distinguishability setting, not a literal percentage of information measured in an experiment.

The app does not demonstrate a derivation of physical matter or spacetime, experimentally confirm GWFM/SEFI/DEFI, or replace quantum mechanics. See [the physics notes](docs/PHYSICS.md) for the formula, assumptions, sampling, and omissions.

## Repository structure

```text
SEFI-Double-Slit-Explorer/
├── public/favicon.svg
├── src/
│   ├── assets/mark.svg
│   ├── components/AdvancedGeometry.tsx
│   ├── components/Concepts.tsx
│   ├── components/Slider.tsx
│   ├── pages/App.tsx
│   ├── simulations/ExperimentCanvas.tsx
│   ├── simulations/physics.ts
│   ├── simulations/physics.test.ts
│   ├── content.ts
│   ├── main.tsx
│   └── styles.css
├── docs/PHYSICS.md
├── docs/CONTENT_GUIDE.md
├── docs/VALIDATION.md
├── .gitignore
├── .nvmrc
├── CONTRIBUTING.md
├── INSTALL.md
├── LICENSE
├── README.md
├── index.html
├── package.json
├── package-lock.json
├── tsconfig.json
└── vite.config.ts
```

Tailwind v4 is configured through the Vite plugin and the `@theme` block in `src/styles.css`; a legacy Tailwind configuration file is not needed. Content lives in `src/content.ts`; the probability model lives separately in `src/simulations/physics.ts`.

## Research and references

- [Feynman Lectures, Volume III, Chapter 1](https://www.feynmanlectures.caltech.edu/III_01.html): quantum interference and the double-slit experiment.
- [Englert, Physical Review Letters 77, 2154 (1996)](https://doi.org/10.1103/PhysRevLett.77.2154): distinguishability and visibility bound.
- [David Tong's quantum field theory lectures](https://www.damtp.cam.ac.uk/user/tong/qft.htm): established field-theory context.
- [Project manuscripts](https://github.com/JasonSEFIDEFI/PhD): the proposed research interpretations.
- [SEFI-PY](https://github.com/JasonSEFIDEFI/SEFI-PY): separate computational research.

See [CONTRIBUTING.md](CONTRIBUTING.md) before changing scientific claims. The included license permits educational evaluation; it does not declare the owner's research public domain or grant a broad commercial license.
