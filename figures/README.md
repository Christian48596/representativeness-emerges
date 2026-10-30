# Reproduce the revised figures

These scripts change presentation or derive summaries from archived numerical
outputs; they do not alter the three original analysis programs. Run the
commands from the repository root. Install the repository requirements first.

## Figures 1 and 2

Extract the original Zenodo version
[10.5281/zenodo.21717613](https://doi.org/10.5281/zenodo.21717613) so the data
folder contains `BSE/` and `reference_outputs/`, then run:

```bash
python figures/regenerate_figure1_layout.py --data-root /path/to/extracted/data
python figures/regenerate_figure2_layout.py --data-root /path/to/extracted/data
```

Alternatively set `REA_DATA_ROOT` to the extracted data folder. If neither an
explicit argument nor this environment variable is provided, the scripts use
`data/` in the repository root. Output defaults to `outputs/figures/`; override
it with `--output-dir /path/to/output`.

Both figures have four identical plot areas, no overall or panel titles, and
bold panel letters outside the upper-left corners. They preserve the complete
original BSE fields of view. The BSE raster is reduced only for display; it is
not resampled for any numerical calculation. Figure 1 uses image 4, its
archived binary mask, and its window-diagnostic CSV. Figure 2 uses image 1 and
`derived_data/figure2_plot_data.npz`. Its centred square indicates the selected
support size, not necessarily the location of the admitted spectral window.
Figure 2 retains the original logarithmic axes and labels them by the quantities
being plotted; its spectra are individually peak-normalized for display only.

PNG output defaults to 600 dpi. PDF output retains vector text, axes, and curves
with raster image panels. Figure 2 accepts `--dpi` for lower-resolution previews.
Both scripts write a layout-audit JSON with panel dimensions.

## Figure 3 and supplementary property maps

Using the bundled compact inputs only:

```bash
python figures/regenerate_revised_figures.py
```

If `REA_DATA_ROOT` is set, this script instead reads the original archived maps
and cell responses. To request that explicitly:

```bash
python figures/regenerate_revised_figures.py --data-root /path/to/extracted/data
```

Output defaults to `outputs/revised_figures/IMG/` and
`outputs/revised_figures/derived_data/`. The `--output-dir` option changes their
common parent. The script regenerates Figure 3, its individual panels, the
reconstructed QEMSCAN map figure, and derived CSV summaries. Figure 3(a) pools
pixels from all seven complete property maps and uses a linear percentage axis.
The elastic summaries invert each full symmetric 3-by-3 plane-strain stiffness
matrix, retaining normal–shear coupling, before calculating `Ex`, `Ey`, and
`Gxy`. Pooled-window SEMs are descriptive; they are not independent-image
confidence intervals. Open symbols at 384 property pixels denote image 6 only.
