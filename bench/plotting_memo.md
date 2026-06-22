# Plotting memo — figure workflow

Working reference for reworking the figures and producing notebooks. **Project context:** the appendix
builds with **lualatex + OldStandard** fonts; current figures are **gnuplot SVG → `rsvg-convert` → PDF**
from `.dat` tables in `bench/`; the comparison plot uses **matplotlib + seaborn**.

## 1. Notebooks (development + interactive/web)

- **Authoring:** **marimo** (pure `.py`, git-friendly, reactive, no hidden Jupyter state). Notebooks
  orchestrate only — heavy logic stays in `vdjmatch`/`bench` modules; fixed seeds; strip large outputs
  before commit.
- **Docs gallery:** Jupyter `.ipynb` via **nbsphinx** (`nbsphinx_execute=never`; run in VS Code, commit
  outputs) — already wired in `docs/`.
- **Interactive / web branch:** marimo **apps** — `mo.ui` sliders + `mo.hstack/vstack/grid` → dashboards,
  exported to **WASM/HTML** for a server-less web demo. Interactive charts: **altair** or **plotly**
  (render in marimo/Jupyter, export to HTML).
- Every notebook ends by **saving the static vector figure** (cases 2/3) so paper assets are
  reproducible from the same code. `pip install marimo altair`

## 2. Standalone vector figure from a table → LaTeX

One data table (`bench/*.dat|tsv`) → one vector figure; figure ← table ← bench script (commit both).

- **Renderer:** keep **gnuplot** for simple line/bar (dependency-light, SVG/PDF terminals); **matplotlib**
  for anything richer.
- **LaTeX-perfect text** (match the OldStandard doc font, not Helvetica): render via the matplotlib
  **`pgf` backend** (text typeset by lualatex with the doc preamble) or gnuplot's `cairolatex`/`epslatex`
  terminal. Robust fallback stays **SVG → `rsvg-convert` → PDF → `\includegraphics`** (current).
- Pin a tiny shared style (font size, line width, palette) so all standalone figures match.

```python
import matplotlib as mpl
mpl.use("pgf"); mpl.rcParams.update({"pgf.texsystem": "lualatex",
    "pgf.preamble": r"\usepackage{OldStandard}\usepackage{unicode-math}"})
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(3.3, 2.4))    # final column width
# ... plot from the table ...
fig.savefig("appendix/figures/foo.pdf")        # text now matches the document
```

## 3. Assembled panels + harmonization → LaTeX

Multi-panel composites (manuscript Fig 1–4): several sub-panels, consistent fonts/colors/sizes, panel
labels A/B/C, drawn at final journal size.

- **Layout:** matplotlib **`subplot_mosaic`/GridSpec + `constrained_layout=True`** (no extra dep), or
  **UltraPlot** (auto tight layout, shared axes, colorbars/legends) for complex grids; **seaborn** for
  panel content. Compose separately-generated panels with **svgutils** (consistent A/B/C labels);
  **Veusz** is the GUI option for visual assembly.
- **Harmonization:** ONE global style file + ONE palette (`seaborn.color_palette("colorblind")`),
  perceptually-uniform colormaps (`cmasher`/`colorcet`); set everything at final size.
- **LaTeX:** export the whole composite as **one vector PDF** → `\includegraphics` (internal A/B/C
  labels), or the `pgf` backend for text-perfect fonts. Use `subcaption`/`\subfloat` only when panels are
  separate files; prefer one composed figure. `pip install matplotlib seaborn ultraplot svgutils cmasher`

## Defaults for this project
- **Notebooks:** marimo (+ nbsphinx for docs); interactive demos exported to HTML/WASM.
- **Standalone appendix figures:** gnuplot now → **matplotlib `pgf`** for the final OldStandard-matched set.
- **Manuscript composites:** matplotlib/UltraPlot + seaborn, one style file, one vector PDF per figure.
