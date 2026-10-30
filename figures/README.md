# Manuscript figures

The three scripts in this directory assemble the complete figures used in the
manuscript from the analysis outputs. They set the panel arrangement, lettering,
axes and legends and export both PNG and PDF files. They do not estimate a new
representative support or solve a cell problem.

| Script | Assembled figure | Required data |
| --- | --- | --- |
| `figure_1.py` | `Figure_1.png` and `Figure_1.pdf`: image, admissible domain and screening diagnostics | Image, matching admissibility mask and detector window CSV. |
| `figure_2.py` | `Figure_2.png` and `Figure_2.pdf`: structural convergence and covariance spectra | Representative image and plotting NPZ. |
| `figure_3.py` | `Figure_3.png` and `Figure_3.pdf`: phase composition and apparent physical properties | Property maps and per-window solver responses, or the bundled study inputs. |

All commands below run from the repository root after installing
`requirements.txt`. Each script accepts `--help`, `--output-dir` and `--dpi`.
The default output directory for all three scripts is `outputs/figures/`.
Use a separate output directory for each dataset.

## Generate Figures 1, 2 and 3 for the manuscript

Download and extract the study data as described in the
[repository README](../README.md). The data directory contains `BSE/`,
`QEMSCAN/` and `reference_outputs/`. Then run:

```bash
export REA_DATA_ROOT="/absolute/path/to/extracted/data"

python figures/figure_1.py --data-root "$REA_DATA_ROOT"
python figures/figure_2.py --data-root "$REA_DATA_ROOT"
python figures/figure_3.py --data-root "$REA_DATA_ROOT" --study-preset
```

Figure 1 uses BSE image 4 and its screening diagnostics. Figure 2 uses BSE image
1 and the bundled `derived_data/figure2_plot_data.npz` curves.
`--study-preset` selects the study's mineral-property assignments, 0.01 mm
property-map calibration and nominal 2.04 mm support marker for Figure 3.

Figure 3 and its property-map illustration can also be generated entirely from
the compact inputs included in the repository:

```bash
env -u REA_DATA_ROOT python figures/figure_3.py
```

This command uses the bundled phase counts, image-1 property maps and per-window
cell responses. Figures 1 and 2 require the BSE image files. To use newly
calculated structural curves, pass their NPZ explicitly with `--data-file`.
Changing `--data-root` alone does not change Figure 2's numerical curve file.

## Figure 1: image and admissibility panels

For a field with any filename:

```bash
python figures/figure_1.py \
  --image-file data/my_field.tif \
  --mask-file outputs/my_masks/my_field_stationary_mask.tif \
  --windows-file outputs/my_masks/my_field_stationary_windows.csv \
  --output-dir outputs/my_dataset/figures
```

The mask must match the image's pixel dimensions. The detector CSV must contain
a complete rectangular grid with columns `iy`, `ix`, `mean`, `mu_residual`,
`std_residual` and `stationary_clean`. The domain detector produces these columns.
The four panels have equal dimensions and bold external letters, with no panel
titles. Ticks adapt to the image and grid dimensions. The displayed raster is
reduced only for rendering; the analysis data are not resampled by the plotter.

Without explicit input paths, the script uses `BSE/4_cut.tif` and its matching
mask and CSV under
`reference_outputs/stationary_domains_FULLRES_R2048_STRICT/`.
The data root is selected by `--data-root`, then `REA_DATA_ROOT`, then the
repository's `data/` directory. Explicitly supply all three input paths for a
custom dataset.

## Figure 2: structural convergence panels

```bash
python figures/figure_2.py \
  --image-file data/my_field.tif \
  --data-file outputs/my_validation/figure2_plot_data.npz \
  --output-dir outputs/my_dataset/figures
```

The image is a representative field from the analysed ensemble. The NPZ is
created by `validation/build_plot_data.py` and can describe any number of images
and any tested integer support sequence. It contains `L`, `sample_mean_Z`,
`eta_C`, `eta_C_sem`, `tau_Z`, `tau_C`, `accepted_L` and `reference_L`, plus
optional `k_L`, `p_L` and `sem_L` spectrum arrays for each displayed support.
NPZ files are loaded with `allow_pickle=False`.

When no support is selected, `accepted_L` is numeric NaN and the figure contains
no selection marker or image square. If no common admissible supports exist,
the panels state that explicitly. Spectrum keys are discovered from the file.
The four panels have equal dimensions and bold external letters.

Support and wavenumber axes use pixels and inverse pixels by default. Supply
`--pixel-size-mm` or include a positive `pixel_size_mm` scalar in the NPZ to
convert those axes to millimetres and inverse millimetres. This conversion
assumes the analysis uses `dx=1` in pixel coordinates. Image-coordinate axes
remain in pixels. The white square illustrates the selected size; it does not
identify the sampled spectral-window location.

Spectra are peak-normalized for display only; the structural selection uses
unnormalized spectra. Without explicit paths the script uses the study's
`BSE/1_cut.tif` and `derived_data/figure2_plot_data.npz`.

## Figure 3: phase composition and physical response

For any number of matched property maps with standard filenames:

```bash
python figures/figure_3.py \
  --data-root /absolute/path/to/my_dataset \
  --window-values outputs/my_property_run/window_values.csv \
  --property-pixel-size-mm 0.01 \
  --rea-pixels 160 \
  --output-dir outputs/my_dataset/figures
```

Replace the example calibration and support with the values for your dataset.
The data root contains `QEMSCAN/kappaMap_ID.csv`, `EMap_ID.csv` and
`nuMap_ID.csv`, or corresponding `.npy` files. Matching suffixes identify a map
triplet. The script discovers all matching sets, including nine, twenty or more,
and orders the suffixes naturally. Each triplet must have matching dimensions.

For other filenames, supply explicit lists in matching order:

```bash
python figures/figure_3.py \
  --kappa-maps data/sampleA_kappa.csv data/sampleB_kappa.csv \
  --E-maps data/sampleA_E.csv data/sampleB_E.csv \
  --nu-maps data/sampleA_nu.csv data/sampleB_nu.csv \
  --window-values outputs/my_property_run/window_values.csv \
  --property-pixel-size-mm 0.01 \
  --output-dir outputs/my_dataset/figures
```

The three lists must have equal lengths. `--map-index` chooses the one-based
list position used in the map illustration. Supply `--rea-mm` or `--rea-pixels`
when a selected support is available; the latter requires
`--property-pixel-size-mm`. Omit these markers when no support is selected.

`--window-values` accepts the solver's per-window CSV. Required columns are
`image`, `window_index`, `L`, `K11`, `K22`, `C11`, `C12`, `C13`, `C22`, `C23` and
`C33`. Supply `L_mm` in the table or a pixel calibration. The table may contain
only the subset of images with accepted windows, and image identifiers may be
strings. The script inverts the full stiffness matrix, retaining coupling, to
calculate `Ex`, `Ey` and `Gxy`. Existing `Ex` and `Ey` columns are checked against
that inversion.

Open symbols identify supports with a reduced response-image population. The
image sets are determined from the CSV. Single-window SEM is undefined and
written as NaN. Coupling components that include nonpositive values use a
symmetric-log scale.

Only the response CSV is required to generate the response panels:

```bash
python figures/figure_3.py \
  --window-values outputs/my_property_run/window_values.csv \
  --output-dir outputs/my_dataset/figures
```

### Phase identification

Without a phase table, Figure 3 contains three response panels and the map
illustration shows conductivity, Young's modulus and Poisson ratio. Generic
property maps receive no assumed mineral identities.

To include a composition panel and classified map, supply `--phase-table
phases.csv` or a JSON table with `name`, `kappa`, `E` and `nu` fields:

```csv
name,kappa,E,nu
Material A,2.0,40.0,0.25
Material B,5.0,100.0,0.30
```

Use the units of the numerical maps: conductivity in W m^-1 K^-1, Young's modulus
in GPa and dimensionless Poisson ratio. A blank CSV value or JSON null is a
wildcard for that property. Every pixel must match exactly one class; unmatched
or ambiguous assignments produce an error. Composition uses area-weighted pixel
counts over all supplied maps. Classification identifies phases from assigned
property values; it does not reconstruct mineral classification from microscopy
images.

Use `--study-preset` only for maps with the dune-sand study's exact mineral-property
assignments. The argument-free compact-data route selects those assignments
automatically. An explicit data root or `REA_DATA_ROOT` selects a dataset and
requires an explicit preset or phase table for classification.

## Output files and rendering

The assembled `Figure_1`, `Figure_2` and `Figure_3` PNG/PDF files are written
directly to the selected output directory. Figure 3 also exports its individual
panels; for the manuscript these are `Figure_3a`, `Figure_3b`, `Figure_3c` and
`Figure_3d`. Without a phase table, the three response panels are
`Figure_3a`, `Figure_3b` and `Figure_3c`. Its map illustration is `property_maps.png` and
`property_maps.pdf`. Numerical summaries and Figure 3 metadata are stored in the
`derived_data/` subdirectory of the output directory.

When classification is omitted, the script removes phase-summary CSVs and any
inapplicable fourth-panel export from the selected destination. When maps are
omitted, it removes its map
illustration. Use a separate destination for each dataset to keep outputs
associated with their inputs.

Figures 1 and 2 default to 600 dpi; Figure 3 defaults to 400 dpi. Use `--dpi` for
lower-resolution previews. PDF exports preserve vector axes and text and embed
raster image panels.
