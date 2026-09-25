# Installation and sharing

## Development

1. Install Node.js 24 LTS (the version family recorded in `.nvmrc`).
2. Extract or clone the complete project, keeping `package.json` in the root.
3. Open a terminal in that folder.
4. Run `npm install` on the first installation. After a lockfile is present, `npm ci` reproduces its dependency versions.
5. Run `npm run dev` and visit http://127.0.0.1:5173.

The development server binds to this computer only. If port 5173 is occupied, stop the other server or use `npm run dev -- --port 5174` and open that address. A strict port setting prevents quietly opening the wrong app.

## Production preview

```sh
npm test
npm run build
npm run preview
```

Use the address printed by the preview server, normally http://127.0.0.1:4173. Do not open `dist/index.html` directly with a file URL: serve the directory over HTTP.

## Public demo

The build is a static website. Upload the **contents of `dist/`** to a static host. Vite uses a relative asset base so a repository subdirectory can work without an absolute site root. There are no client-side routes requiring rewrite rules.

A local `127.0.0.1` address works only on your own computer; it is not a public share link. Publish the static build before sharing a demo on LinkedIn or with a reviewer. Include the research limitations and source references with any public presentation. No hosting account or repository is created by these instructions.

## Browser requirements

A modern browser supporting Canvas 2D, ResizeObserver, requestAnimationFrame, and JavaScript modules. If Canvas is unavailable, the app displays a fallback explanation. For a lower-motion demonstration, use Pause or enable reduced motion in the operating system before opening the app.

## Isolated from existing software

This application has its own dependencies and launch command. Do not install these JavaScript dependencies inside the SEFI-PY engine environment. No Python engine configuration or source file is needed or modified.
