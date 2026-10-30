# Figure generation for study and custom datasets

These scripts plot supplied numerical data. They do not estimate a new REA or
solve a cell problem. Run the analysis first, then select the corresponding
image, diagnostic CSV, plotting NPZ or per-window property CSV.

All commands below run from the repository root. The required packages are in
`requirements.txt`. Each script accepts `--help`, `--output-dir` and `--dpi`.
Use a separate output directory for each dataset.

## Image and admissibility panels: Figure 1

For an image with any filename:

```bash
python figures/regenerate_figure1_layout.py \
  --image-file data/my_field.tif \
  --mask-file outputs/my_masks/my_field_stationary_mask.tif \
  --windows-file outputs/my_masks/my_field_stationary_windows.csv \
  --output-dir outputs/my_dataset/figures
```

The mask must match the image's pixel dimensions. The detector CSV must contain
a complete rectangular grid with columns `iy`, `ix`, `mean`, `mu_residual`,
`std_residual` and `stationary_clean`. These are the columns produced by the
core detector. The four panels have equal dimensions and bold external letters;
there are no panel titles. Ticks adapt to the image and grid sizes. The displayed
raster is reduced only for rendering; the plotter does not alter analysis data.

With no explicit input paths, the plotter uses `BSE/4_cut.tif` under its data
root and the corresponding diagnostics in
`reference_outputs/stationary_domains_FULLRES_R2048_STRICT/`. The data root is
selected by `--data-root`, then `REA_DATA_ROOT`, then the repository's `data/`.
For a custom dataset, explicitly supplying all three input paths avoids relying
on study-specific directory names.

## Structural convergence panels: Figure 2

```bash
python figures/regenerate_figure2_layout.py \
  --image-file data/my_field.tif \
  --data-file outputs/my_validation/figure2_plot_data.npz \
  --output-dir outputs/my_dataset/figures
```

The image is one representative field from the analysed ensemble. The NPZ is
created by `validation/rebuild_derived.py` and can describe any number of images
and any tested integer support sequence. It contains `L`, `sample_mean_Z`,
`eta_C`, `eta_C_sem`, `tau_Z`, `tau_C`, `accepted_L` and `reference_L`, plus
optional `k_L`, `p_L`, `sem_L` spectrum arrays for each displayed support.

`accepted_L` is numeric NaN when no support was selected. Such data are plotted
without a selection marker or image square. If no common admissible supports
exist, the panels state that explicitly. Spectrum keys are discovered from the
file rather than restricted to the study's four example sizes. NPZ arrays are
loaded with `allow_pickle=False`.

The support and wavenumber axes use pixels and inverse pixels by default. Supply
`--pixel-size-mm` or include a positive `pixel_size_mm` scalar in the NPZ to
convert those axes to millimetres and inverse millimetres. This conversion
assumes that the input wavenumbers use the core workflow's `dx=1` pixel
convention. Image-coordinate axes remain in pixels. The white square is a
schematic of the selected size, not the admitted spectral-window location.

The spectra are peak-normalized for display only; the structural analysis uses
unnormalized spectra. The plotter preserves the input numerical curves. Without
explicit paths it uses the study's `BSE/1_cut.tif` and bundled
`derived_data/figure2_plot_data.npz`; changing `--data-root` alone does not change
the bundled numerical curve file, so use `--data-file` for new calculations.

## Property responses and maps: Figure 3 and map illustration

For an arbitrary number of matched property maps with standard filenames:

```bash
python figures/regenerate_revised_figures.py \
  --data-root /absolute/path/to/my_dataset \
  --window-values outputs/my_property_run/window_values.csv \
  --property-pixel-size-mm 0.01 \
  --rea-pixels 160 \
  --output-dir outputs/my_dataset/property_figures
```

The example calibration and selected size must be replaced with the values for
your dataset. The data root contains `QEMSCAN/kappaMap_ID.csv`, `EMap_ID.csv`
and `nuMap_ID.csv` (or `.npy` files), where matching suffixes identify a map set.
The script discovers every matching set, including nine, twenty or more sets.
It requires matching suffixes across the three fields and orders them naturally.

For other filenames, supply explicit lists in matching order:

```bash
python figures/regenerate_revised_figures.py \
  --kappa-maps data/sampleA_kappa.csv data/sampleB_kappa.csv \
  --E-maps data/sampleA_E.csv data/sampleB_E.csv \
  --nu-maps data/sampleA_nu.csv data/sampleB_nu.csv \
  --window-values outputs/my_property_run/window_values.csv \
  --property-pixel-size-mm 0.01 \
  --output-dir outputs/my_dataset/property_figures
```

Each map triplet must have matching shapes. The three lists must have the same
length; any number of triplets may be supplied. `--map-index` chooses the
one-based list position shown in the map illustration. Omit an REA marker when
no support was selected; use `--rea-mm` or `--rea-pixels` when a selected support
is available. The latter requires `--property-pixel-size-mm`.

`--window-values` accepts the core solver's per-window CSV. Required columns are
`image`, `window_index`, `L`, `K11`, `K22`, `C11`, `C12`, `C13`, `C22`, `C23`, and
`C33`. Supply `L_mm` or a pixel calibration. The response table may contain only
the subset of images with accepted windows, and image identifiers may be strings.
The script inverts each full stiffness matrix, retaining coupling, to calculate
`Ex`, `Ey` and `Gxy`. If existing `Ex` and `Ey` columns are present, they are
checked against that inversion. No stiffness components are discarded.

The plots distinguish supports with the complete response-image population
from those with a reduced image set. These sets are determined from the CSV;
no image count or image identifier is fixed. Single-window SEM is undefined
and written as NaN. Negative coupling components are shown on a symmetric-log
scale instead of silently disappearing from a logarithmic plot.

Only the response CSV is required for the response plots:

```bash
python figures/regenerate_revised_figures.py \
  --window-values outputs/my_property_run/window_values.csv \
  --output-dir outputs/my_dataset/property_figures
```

### Optional phase identification

Generic property maps are not automatically assigned dune-sand mineral names.
Without a phase table, Figure 3 contains the three response panels and the map
illustration shows conductivity, Young's modulus and Poisson ratio.

To add a composition panel and classified map, supply `--phase-table phases.csv`
or a JSON table. The table has `name`, `kappa`, `E` and `nu` fields:

```csv
name,kappa,E,nu
Material A,2.0,40.0,0.25
Material B,5.0,100.0,0.30
```

Use physical values consistent with the input maps: conductivity in W m^-1 K^-1,
Young's modulus in GPa and dimensionless Poisson ratio. A blank CSV value or
JSON null is a wildcard for that property. Every pixel must match exactly one
class; unmatched or ambiguous classifications produce an error. Composition
uses area-weighted pixel counts over all supplied maps. This optional step
identifies classes from already assigned property values; it does not recover
upstream mineral classification from microscopy images.

The plotter removes only its own stale phase-chart files when phase
classification is omitted and its own stale map illustration when no maps are
supplied, preventing previous classifications from appearing as current output.

## Reproduce the dune-sand figures

Use the extracted published data directory for Figures 1 and 2:

```bash
python figures/regenerate_figure1_layout.py --data-root /path/to/extracted/data
python figures/regenerate_figure2_layout.py --data-root /path/to/extracted/data
python figures/regenerate_revised_figures.py \
  --data-root /path/to/extracted/data --study-preset
```

`--study-preset` explicitly selects the published mineral assignments, the
0.01 mm property-pixel calibration and the nominal 2.04 mm support marker.
The compact Figure 3/SI route needs no downloaded maps:

```bash
env -u REA_DATA_ROOT python figures/regenerate_revised_figures.py
```

This argument-free data route uses the bundled study maps and cell responses.
The numerical derived CSVs retain the study's reference values. `--data-root`
or `REA_DATA_ROOT` otherwise selects a dataset; supply `--study-preset` only if
its maps use the study's exact mineral-property assignments.

## Outputs

Figures 1 and 2 default to `outputs/figures/`. Property figures default to
`outputs/revised_figures/IMG/`, with derived CSVs and provenance metadata under
`outputs/revised_figures/derived_data/`. Existing `_revised` filenames are
retained for compatibility. Figures 1 and 2 default to 600 dpi; property figures
default to 400 dpi. Each script can render lower-resolution previews using
`--dpi`. PDF files preserve vector axes/text and embed the raster image panels.
