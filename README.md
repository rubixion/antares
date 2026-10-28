# Antares

**Live:** https://rubixion.github.io/antares/

An interactive 3D sky in the browser. It shows all 88 constellations and every naked-eye star. Click a star to fly into its planetary system, which is built from every confirmed exoplanet in the NASA Exoplanet Archive.

## Features

- Every star down to magnitude 6.5, with constellation figures, IAU boundaries, labels, the Milky Way, and an RA/Dec grid
- Click a star for its details, then click again to fly into its system
- Planetary systems use real orbital elements where they're known
- Multiple-star systems come from Eggleton & Tokovinin (2008), with measured binary orbits from ORB6 and SB9
- The Solar System's positions come from JPL Horizons, as of the date the data was built
- Optional simulated planets for stars with no confirmed ones (they're clearly labelled as simulated)
- Search for a constellation (*Orion*), a star (*Betelgeuse*, *HIP 27989*), a planet host (*TRAPPIST-1*), or a planet (*Kepler-22 b*)
- Settings for each sky layer, orbit speed, auto-rotate, and the background

## Controls

| Input | Action |
|---|---|
| Drag | Look around / orbit |
| Scroll | Zoom |
| Click star | Star info, and click again to fly in |
| Click sky | Select that constellation |
| Click planet | Planet info, and follow it |
| Enter | Fly into the selected system |
| Esc | Close, or go back to the sky |
| R | Reset view |
| / | Search |

## Run locally

There's no build step. Serve the folder with any static server:

```sh
npx serve .
# or
python -m http.server
```

Opening `index.html` directly as a `file://` page won't work, because browsers block ES module imports from file URLs.

## Rebuilding the data

`data.js` is generated, so don't edit it by hand. To refresh the catalog snapshot:

```sh
pip install numpy
python tools/build_data.py
```

Downloads are cached in `tools/.cache/` (gitignored). Delete that folder to pull fresh data.

## Deploy

It's hosted on GitHub Pages from the root of the `master` branch (Settings → Pages → Deploy from a branch → `master` / `(root)`).

## Stack

- Three.js r169: the WebGPU renderer, loaded from the jsDelivr CDN
- A single `index.html` plus the generated `data.js`
- A Python and NumPy build script for the data

## Data sources

| Data | Source | License |
|---|---|---|
| Stars | [HYG v4.1](https://github.com/astronexus/HYG-Database) | CC BY-SA 4.0 |
| Figures, names, IAU boundaries | [Stellarium](https://github.com/Stellarium/stellarium) | GPL-2.0 |
| Exoplanets | [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu) (pscomppars) | Public |
| Star multiplicity | [Eggleton & Tokovinin 2008](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/MNRAS/389/869) | VizieR |
| Visual binary orbits | [USNO ORB6](https://crf.usno.navy.mil/data_products/WDS/orb6/orb6orbits.txt) | Public |
| Spectroscopic binary orbits | [SB9](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=B/sb9) | VizieR |
| Solar System positions | [JPL Horizons](https://ssd.jpl.nasa.gov/horizons/) | Public |

Values marked `~` in the app are estimates.
