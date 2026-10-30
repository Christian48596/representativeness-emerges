# Statistical admissibility and long-wavelength structural convergence determine representative support in granular materials

**Christian Tantardini, Fernando Alonso-Marroquín, Abdullah Alqubalee and Eduardo Garzanti**

This repository provides a Python workflow for determining the spatial support over which the structure of a granular material can be represented consistently and for testing that support against independent physical properties. It combines local statistical screening, convergence of image means and long-wavelength spatial correlations, and periodic thermal and elastic homogenization.

The accompanying study applies this framework to two-dimensional polished sections of dune sand. Backscattered-electron (BSE) images provide the structural information. Numerical property maps derived from QEMSCAN mineralogical classifications provide independent conductivity and elasticity fields. The repository contains the analysis programs, reference statistics, uncertainty and sensitivity calculations, and scripts for producing the manuscript figures.

The workflow accepts **N two-dimensional images**, including 9, 20 or more, without editing Python source. Input lists determine N. The uncertainty, sensitivity and plotting programs use the supplied dataset and its parameters. The bundled seven-image dune-sand results provide a reproducibility reference.

## Contents

- [Scientific workflow](#scientific-workflow)
- [Repository organization](#repository-organization)
- [Installation and verification](#installation-and-verification)
- [Reproduce the dune-sand study](#reproduce-the-dune-sand-study)
- [Analyse N new images](#analyse-n-new-images)
- [Uncertainty and sensitivity for N images](#uncertainty-and-sensitivity-for-n-images)
- [Independent property calculations](#independent-property-calculations)
- [Generate figures for your dataset](#generate-figures-for-your-dataset)
- [Read and interpret the results](#read-and-interpret-the-results)
- [Troubleshooting](#troubleshooting)
- [Reproducibility, citation and license](#reproducibility-citation-and-license)

## Scientific workflow

A representative support must be associated with a defined material domain, measured field and convergence criterion. The workflow separates three questions.

1. **Statistical admissibility:** where can comparable local statistics be sampled? A moving-window analysis compares local means and standard deviations with robust image-level references. Morphological processing produces a binary admissibility mask.
2. **Structural convergence:** at what support do the sampled mean and the low-wavenumber covariance spectrum stabilize relative to a finite reference? Candidate square windows are evaluated within the admissible domain, and the first support satisfying the joint persistent-convergence rule is selected.
3. **Independent physical assessment:** do conductivity and elastic responses at the transferred support show compatible scale dependence? The structural support is converted to the property-map grid, where periodic finite-element problems evaluate the apparent material response.

These stages distinguish spatial organization from the convergence of a single average. The BSE gray level is an image-contrast field; physical conductivity and elasticity are supplied separately through numerical maps.

The implemented calculations operate on **two-dimensional fields**. The representative elementary area (REA) is reported through the **side length** of its square support, `L_REA`; its area is `L_REA²`. Thus a reported value in millimetres denotes a side length. The corresponding three-dimensional concept is a representative elementary volume (REV), but this repository does not contain a validated three-dimensional solver. The study's mechanical calculations use **plane strain**, and its conduction calculations are in-plane.

## Repository organization

| File or directory | Role |
| --- | --- |
| [`detect_stationary_domains_general.py`](detect_stationary_domains_general.py) | Detect statistically admissible domains in image fields or combined numerical property maps. |
| [`rea_stationary_mask.py`](rea_stationary_mask.py) | Evaluate structural statistics and spectra across support sizes and select a persistent structural support. |
| [`property_full_homogenized_inside_property_masks.py`](property_full_homogenized_inside_property_masks.py) | Solve periodic thermal and elastic cell problems inside admissible property-map domains. |
| [`validation/`](validation/) | Bootstrap, sensitivity, numerical reconstruction and verification scripts. |
| [`validation/checkpoints/`](validation/checkpoints/) | Compressed per-image statistics and spectra for reproducing the uncertainty analysis without processing the large images again. |
| [`validation/results/`](validation/results/) | Reference bootstrap, sensitivity and per-image results. |
| [`validation/inputs/`](validation/inputs/) | Archived structural summary used for numerical comparison. |
| [`figures/figure_1.py`](figures/figure_1.py) | Assemble manuscript Figure 1: image, admissibility mask and screening diagnostics. |
| [`figures/figure_2.py`](figures/figure_2.py) | Assemble manuscript Figure 2: structural mean, convergence residuals and spectra. |
| [`figures/figure_3.py`](figures/figure_3.py) | Assemble manuscript Figure 3: phase fractions and thermal/elastic responses; also export the property-map illustration. |
| [`derived_data/`](derived_data/) | Compact numerical plotting inputs, phase fractions and property summaries. |
| [`docs/STUDY_WORKFLOW.md`](docs/STUDY_WORKFLOW.md) | Complete study workflow, from raw data to assembled manuscript figures. |
| [`tests/test_n_images.py`](tests/test_n_images.py) | Synthetic 9- and 20-image integration checks and input/error cases. |
| [`requirements.txt`](requirements.txt) | Required Python packages. |
| [`requirements-pinned.txt`](requirements-pinned.txt) | Pinned package versions used to verify the numerical reconstruction. |

Large images, full-resolution masks and numerical property maps are distributed through Zenodo. Generated outputs are written under `outputs/` by the commands below.

## Installation and verification

The programs require Python, NumPy, SciPy, Matplotlib and Pillow. The recorded verification environment uses Python 3.12.14. Commands below use **Bash on Linux or macOS**; on Windows, use WSL or adapt the shell syntax. Run commands from the repository root with the virtual environment active.

```bash
git clone https://github.com/Christian48596/representativeness-emerges.git
cd representativeness-emerges
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

For the recorded package versions, use `python -m pip install -r requirements-pinned.txt` in a Python 3.12 environment.

Check the bundled scientific reference and the N-image interfaces:

```bash
python validation/check_reference_results.py
python -m unittest discover -s tests -p 'test_n_images.py' -v
```

The reference check reconstructs the stored baseline, all 10,000 bootstrap draws, leave-one-image-out results, sensitivity results and Figure 2 arrays. It compares scientific values with the bundled references and independently verifies the three computational source files. It needs no downloaded images. The tests exercise distinct synthetic datasets with 9 and 20 images, matched property maps, thermal and plane-strain elastic solutions, custom figures, incomplete inputs and non-selection cases. Synthetic tests verify software behavior; the physical validation remains the dune-sand study.

Every program accepts `--help`. For example:

```bash
python rea_stationary_mask.py --help
python validation/structural_uncertainty.py --help
python figures/figure_3.py --help
```

N is not fixed by the software, but image size, memory and FEM cost limit practical runs. The structural program holds the image ensemble in memory; the uncertainty and sensitivity programs process one image at a time. Each study BSE field at full resolution contains approximately 205 million pixels, so start with the compact verification on a memory-limited machine.

## Reproduce the dune-sand study

The dataset is available through [Zenodo](https://doi.org/10.5281/zenodo.21717612). Select one record version and download:

| File | Contents |
| --- | --- |
| `BSE_full_resolution_images.tar.gz` | Seven full-resolution BSE images. |
| `QEMSCAN_numerical_property_maps.tar.gz` | Numerical conductivity, Young's-modulus and Poisson-ratio maps. |
| `reference_outputs_and_masks.tar.gz` | Masks, numerical outputs and execution logs. |
| `README_DATA.md` | Dataset documentation. |
| `SHA256SUMS.txt` | Checksums for that version. |

Put the files in `data/`. From that directory verify them with `sha256sum -c SHA256SUMS.txt` on Linux or `shasum -a 256 -c SHA256SUMS.txt` on macOS. Then, from the repository root:

```bash
tar -xzf data/BSE_full_resolution_images.tar.gz -C data
tar -xzf data/QEMSCAN_numerical_property_maps.tar.gz -C data
tar -xzf data/reference_outputs_and_masks.tar.gz -C data
export REA_DATA_ROOT="$PWD/data"
```

The extracted directory contains `BSE/1_cut.tif` through `7_cut.tif`, `QEMSCAN/kappaMap_1.csv` through `kappaMap_7.csv`, corresponding `EMap_*.csv` and `nuMap_*.csv`, and `reference_outputs/`.

Follow **[docs/STUDY_WORKFLOW.md](docs/STUDY_WORKFLOW.md)** for the complete calculation and manuscript-figure assembly with the study settings. That guide includes the complete commands for BSE screening, structural selection, property screening, periodic homogenization and figure assembly. Defaults alone do not reproduce those settings.

To reconstruct the stored uncertainty and sensitivity outputs without the large images:

```bash
python validation/structural_uncertainty.py --aggregate --bootstrap 10000 --verify-reference
python validation/screening_sensitivity.py --aggregate
python validation/build_plot_data.py
```

These commands read the bundled `validation/checkpoints/` and write to `outputs/validation/`. `--verify-reference` is an explicit study-only check; do not use it with a new dataset. `check_reference_results.py` additionally compares all stored reference values.

To reproduce the figures using downloaded study data:

```bash
python figures/figure_1.py --data-root "$REA_DATA_ROOT"
python figures/figure_2.py --data-root "$REA_DATA_ROOT"
python figures/figure_3.py --data-root "$REA_DATA_ROOT" --study-preset
```

`--study-preset` selects the study's mineral assignments, 0.01 mm property calibration and nominal 2.04 mm support marker. Figure 2 uses the bundled numerical curves unless `--data-file` is supplied. A compact property-figure route requires no downloaded maps:

```bash
env -u REA_DATA_ROOT python figures/figure_3.py
```

## Analyse N new images

### Inputs and calibration

Use one 2D grayscale field per image file, with unique filename stems such as `sample_A.tif`, `sample_B.tif`, etc. The documented structural pipeline uses an **8-bit grayscale (`L`) conversion**, raw contrast and full image resolution. Higher-bit-depth or color information is not preserved by that conversion; prepare a justified, consistent grayscale representation before using this workflow. This is not an automatic segmentation or mineral-classification program.

All images pooled into one analysis should represent the intended sampling population and share the same pixel calibration and intensity interpretation. Image dimensions may differ, but the support sequence must be feasible across the ensemble. Overlapping crops are not independent realizations for whole-image bootstrap uncertainty. Keep one input dataset and one parameter configuration in each run directory.

The following example uses the study's numerical settings to show a complete workflow. **Choose the support sizes, screening scales, reference, tolerances and pixel calibration for your own data.** The example requires images and admissible regions large enough to contain the tested supports; copying it onto small images will not produce a meaningful estimate. Support sizes and screening strides are in pixels.

### 1. Select all images and settings

Create a directory containing only the images to analyse. A glob automatically supplies all 9, 20 or N matching files; filenames need not contain numbers. Keep this Bash session open for the subsequent commands.

```bash
INPUT_IMAGES="/absolute/path/to/my_images"
RUN="$PWD/outputs/my_dataset"
mkdir -p "$RUN"
shopt -s nullglob
images=("$INPUT_IMAGES/"*.tif)
if (( ${#images[@]} == 0 )); then
  echo "No TIFF images found; check INPUT_IMAGES and the file extension." >&2
  exit 1
fi
printf 'Using %s images\n' "${#images[@]}"

WINDOWS=256,384,512,640,768,896,1024,1280,1536,1792,2048
REFERENCE=2048
SPECTRA=256,512,1024,2048
RSTAT=2048
STRIDE=512
PIXEL_MM=0.00131
N_SAMPLES=121
TAU_MU=0.02
TAU_SIGMA=0.06
TAU_Z=0.025
TAU_C=0.30
MASK_FRACTION=0.98
```

Replace `*.tif` with the appropriate extension, or supply an explicit ordered array such as `images=("/path/field A.tif" "/path/field B.tif")`. Quoted arrays support spaces in paths. The same order is used for one-based image IDs throughout the run. These input files need not be called `1_cut.tif`, etc.

### 2. Screen admissible domains

```bash
python detect_stationary_domains_general.py "${images[@]}" \
  --input-type image --field-transform raw --max-side 0 \
  --output-dir "$RUN/masks" \
  --stationarity-window "$RSTAT" --stride "$STRIDE" \
  --tau-mu "$TAU_MU" --tau-sigma "$TAU_SIGMA" \
  --min-component-windows 4 --dpi 300
```

For `sample_A.tif`, this produces `sample_A_stationary_mask.tif` and `sample_A_stationary_windows.csv` in `$RUN/masks`, plus other diagnostics. Masks match their image dimensions. Review the retained domain and diagnostic maps before interpreting a structural estimate. The morphological operations remain enabled.

### 3. Determine structural support

```bash
python rea_stationary_mask.py "${images[@]}" \
  --output-dir "$RUN/structural" --output-prefix REA \
  --representative-image 1 --max-side 0 --pixel-size 1 --field-mode raw \
  --stationary-mask-dir "$RUN/masks" \
  --stationary-mask-suffix _stationary_mask.tif \
  --min-mask-fraction "$MASK_FRACTION" \
  --window-sizes "$WINDOWS" --reference-L "$REFERENCE" \
  --selected-spectra "$SPECTRA" --n-common-k-grid 600 \
  --n-window-samples "$N_SAMPLES" --errorbar-mode sem \
  --spectral-band low --low-k-fraction 0.25 --low-k-weight inverse \
  --tau-Z "$TAU_Z" --tau-C "$TAU_C" --tau-ens 0.10 --tail-fraction 0.5 \
  --panel-b-mode raw+smooth --panel-b-smooth-window 3 \
  --panel-c-mode raw+smooth --panel-c-smooth-window 3 --dpi 300 --no-pdf
```

The principal output is `$RUN/structural/REA_combined_window_statistics.csv`. Use a basename for `--output-prefix`. The spectrum remains in pixel coordinates (`--pixel-size 1`); multiply selected side lengths by `PIXEL_MM` for millimetres. The actual reference can fall back to the largest common admissible support if the requested reference is unavailable; inspect the reported reference and available sizes.

The structural program requires enough common admissible supports. If there are too few, it can stop with an informative error. The validation workflow below additionally records insufficient support and non-selection in machine-readable output. An absent selection is a result to report, not a size to replace with the largest window.

## Uncertainty and sensitivity for N images

The uncertainty and sensitivity programs import the core numerical routines. They support the full-resolution, raw grayscale, low-wavenumber protocol shown above; they do not mirror every alternative preprocessing mode of the core command-line program.

Define the analysis settings once and pass the same configuration to both image-processing commands. This example tests screening sides 1536 and 3072, shorter references 1536 and 1792, and larger-reference availability through 4096. Adjust those sizes together with your baseline grid. At least four distinct structural support sizes are required, and the requested reference must be in that grid.

```bash
analysis_options=(
  --images "${images[@]}"
  --mask-dir "$RUN/masks"
  --window-sizes "$WINDOWS" --reference-L "$REFERENCE"
  --selected-spectra "$SPECTRA"
  --n-window-samples "$N_SAMPLES" --n-common-k-grid 600
  --tau-Z "$TAU_Z" --tau-C "$TAU_C"
  --min-mask-fraction "$MASK_FRACTION" --tail-fraction 0.5
  --low-k-fraction 0.25 --low-k-weight inverse
  --extra-reference-sizes 2560,3072,3584,4096
  --short-reference-sizes 1536,1792
  --screening-sizes 1536,3072 --screening-stride "$STRIDE"
  --screening-tau-mu "$TAU_MU" --screening-tau-sigma "$TAU_SIGMA"
  --min-component-windows 4 --pixel-size-mm "$PIXEL_MM"
)

python validation/structural_uncertainty.py --all-images \
  "${analysis_options[@]}" \
  --checkpoint-dir "$RUN/checkpoints" --output-dir "$RUN/validation"

python validation/screening_sensitivity.py --all-images \
  "${analysis_options[@]}" \
  --checkpoint-dir "$RUN/checkpoints" --output-dir "$RUN/validation"

python validation/structural_uncertainty.py --aggregate \
  --checkpoint-dir "$RUN/checkpoints" --output-dir "$RUN/validation" \
  --expected-images "${#images[@]}" --bootstrap 10000 --seed 20261029

python validation/screening_sensitivity.py --aggregate \
  --checkpoint-dir "$RUN/checkpoints" --output-dir "$RUN/validation" \
  --expected-images "${#images[@]}"

python validation/build_plot_data.py \
  --checkpoint-dir "$RUN/checkpoints" --output-dir "$RUN/validation" \
  --expected-images "${#images[@]}"
```

`--all-images` processes N images sequentially. To process one position, replace it with `--image 3` while retaining the **complete input list** and settings. An alternative is `--data-root /path/to/dataset`, which discovers images under `BSE/`; explicit `--images` avoids study-specific directory assumptions. Use `--masks PATH...` instead of `--mask-dir` for explicitly paired mask files.

Always pass the same `--checkpoint-dir` to processing, aggregation and reconstruction. Omitting it during aggregation selects the bundled seven-image checkpoints. A new dataset or configuration needs a new checkpoint directory. Stored metadata records ordered image and mask paths, file sizes and modification times, and configuration. The programs reject incomplete or mixed populations and changed inputs at checkpoint-writing time. These metadata checks are not content checksums. Aggregation reads the stored metadata and can run without access to the original image files.

Each bootstrap draw samples **N images with replacement from the N supplied images**. Supports remain linked within each image. Leave-one-image-out analysis creates N exclusions when N > 1. With one image, independent-image uncertainty cannot be estimated; output marks the degenerate bootstrap and omits leave-one-image-out estimates.

| Generated file in `$RUN/validation` | Contents |
| --- | --- |
| `validation_results.json` | Baseline, masked/unmasked sampling comparisons, availability, bootstrap and leave-one-image-out results. |
| `bootstrap_draws.npz` | Draw indices, selected supports and spectral results. |
| `baseline_per_image_1.csv` … `baseline_per_image_N.csv` | Individual-image structural statistics, when a baseline curve can be calculated. |
| `rstat_sensitivity_results.json` | Configured screening-window experiments. |
| `reference_sensitivity_extra.json` | Shorter-reference and explicitly flagged exploratory larger-reference diagnostics. |
| `figure2_plot_data.npz` | Curves, selected spectra, calibration and selected support for Figure 2. |

A non-selection remains explicit, including when every bootstrap draw fails or no common support is available. Successful-only percentiles must be reported together with the failure fraction. Extra-reference diagnostics that bypass the 49-position availability check are marked `enforces_production49_common_size_gate: false` and must not be presented as production-protocol results. See [validation/README.md](validation/README.md) for defaults and file formats.

## Independent property calculations

This stage is optional when only structural support is being studied. To assess thermal and elastic response, supply **numerical**, registered conductivity, Young's-modulus and Poisson-ratio maps. The three maps within each triplet must have identical shapes and corresponding pixels. Use conductivity in W m^-1 K^-1, Young's modulus in GPa and dimensionless Poisson ratio. Do not use intensity-rescaled illustrations as physical inputs.

The number of property-map triplets is inferred from their lists. It need not equal the structural-image count, but its sampling population and correspondence must be justified. The shell loop below pairs standard filenames by suffix, avoiding accidental list-order mismatches. Change the paths and filenames to match your data; explicit ordered arrays also work.

```bash
PROPERTY_DIR="/absolute/path/to/numerical_maps"
kappa_files=("$PROPERTY_DIR/"kappaMap_*.csv)
elastic_files=()
poisson_files=()
if (( ${#kappa_files[@]} == 0 )); then
  echo "No numerical property maps found." >&2
  exit 1
fi
for kappa in "${kappa_files[@]}"; do
  map_id="${kappa##*/kappaMap_}"
  map_id="${map_id%.csv}"
  efile="$PROPERTY_DIR/EMap_${map_id}.csv"
  nufile="$PROPERTY_DIR/nuMap_${map_id}.csv"
  if [[ ! -f "$efile" || ! -f "$nufile" ]]; then
    echo "Incomplete property-map triplet: $map_id" >&2
    exit 1
  fi
  elastic_files+=("$efile")
  poisson_files+=("$nufile")
done

PROPERTY_PIXEL_MM=0.01
PROPERTY_RSTAT=256
PROPERTY_STRIDE=64
PROPERTY_WINDOWS=32,48,64,80,96,112,128,160,192,204,224,256,320,384,448,512
PROPERTY_REA_PIXELS=204
```

The values above illustrate the study's property settings. Set your own screening scales and support sequence. Convert a structural side to the property grid using `L_property_pixels = L_structural_pixels × structural_pixel_mm / property_pixel_mm`; choose and record an integer rounding convention. The study uses a **nominal 204-pixel (2.04 mm)** transfer. That is not the exact nearest-integer conversion of 1536 × 0.00131 / 0.01, which is 201.216 pixels. For a new dataset, use its measured support and calibration rather than copying 204.

```bash
python detect_stationary_domains_general.py \
  --kappa-maps "${kappa_files[@]}" --E-maps "${elastic_files[@]}" \
  --nu-maps "${poisson_files[@]}" \
  --stationarity-source combined --combined-transform log-robust --max-side 0 \
  --output-dir "$RUN/property_masks" \
  --stationarity-window "$PROPERTY_RSTAT" --stride "$PROPERTY_STRIDE" \
  --tau-mu "$TAU_MU" --tau-sigma "$TAU_SIGMA" \
  --min-component-windows 4 --dpi 300

python property_full_homogenized_inside_property_masks.py \
  --kappa-maps "${kappa_files[@]}" --E-maps "${elastic_files[@]}" \
  --nu-maps "${poisson_files[@]}" \
  --mask-dir "$RUN/property_masks" \
  --mask-template 'property_combined_{i}_stationary_mask.tif' \
  --output-dir "$RUN/properties" --output-prefix properties \
  --window-sizes "$PROPERTY_WINDOWS" --L-rea "$PROPERTY_REA_PIXELS" \
  --property-pixel-size-mm "$PROPERTY_PIXEL_MM" \
  --min-mask-fraction "$MASK_FRACTION" --n-window-samples 4 \
  --elastic-mode plane_strain --errorbar-mode sem \
  --rng-seed 12345 --dpi 300 --no-pdf
```

The detector numbers property masks by one-based **list position**, so use exactly the same list order in both commands. The combined screening field log-transforms and robustly standardizes conductivity and modulus; Poisson ratio is robustly standardized without a logarithm. Positive, physically appropriate coefficients are needed for the cell problems.

The solver uses periodic Q4 elements, 2 × 2 Gauss integration, zero-mean constraints and sparse direct solves for two thermal and three elastic loading cases. The per-window responses are saved to `$RUN/properties/properties_window_values.csv`, alongside a support-wise summary and plots. The mechanical assumption here is **2D plane strain**. The figure postprocessor inverts the full 3 × 3 stiffness matrix to obtain `Ex = 1/S11`, `Ey = 1/S22` and `Gxy = 1/S33`, including directional coupling.

## Generate figures for your dataset

The calculation programs write numerical results and diagnostic plots. The three scripts in `figures/` assemble those results into the manuscript layouts. Each script draws the complete multi-panel figure from numerical inputs, so the manuscript figures can be reproduced without manual image editing. The Figure 3 script also exports individual panels and a property-map illustration.

| Assembly script | Manuscript output | Numerical inputs |
| --- | --- | --- |
| `figures/figure_1.py` | `Figure_1.png` and `Figure_1.pdf` | Selected image, its mask and detector window CSV. |
| `figures/figure_2.py` | `Figure_2.png` and `Figure_2.pdf` | Selected image and structural plotting NPZ. |
| `figures/figure_3.py` | `Figure_3.png` and `Figure_3.pdf` | Property-map triplets, per-window responses and phase assignments. |

Use explicit paths so all panels refer to the intended calculation:

```bash
example_image="${images[0]}"
example_stem="${example_image##*/}"
example_stem="${example_stem%.*}"

python figures/figure_1.py \
  --image-file "$example_image" \
  --mask-file "$RUN/masks/${example_stem}_stationary_mask.tif" \
  --windows-file "$RUN/masks/${example_stem}_stationary_windows.csv" \
  --output-dir "$RUN/figures"

python figures/figure_2.py \
  --image-file "$example_image" \
  --data-file "$RUN/validation/figure2_plot_data.npz" \
  --pixel-size-mm "$PIXEL_MM" --output-dir "$RUN/figures"

python figures/figure_3.py \
  --kappa-maps "${kappa_files[@]}" --E-maps "${elastic_files[@]}" \
  --nu-maps "${poisson_files[@]}" \
  --window-values "$RUN/properties/properties_window_values.csv" \
  --property-pixel-size-mm "$PROPERTY_PIXEL_MM" \
  --rea-pixels "$PROPERTY_REA_PIXELS" \
  --output-dir "$RUN/property_figures"
```

Run the third command only after the property calculations. Omit `--rea-pixels` when no structural support was selected. Supplying only `--window-values` produces response plots without requiring maps. To add phase-composition plots, supply an explicit `--phase-table` with `name,kappa,E,nu` entries; see [figures/README.md](figures/README.md). Generic maps receive no automatic dune-sand mineral names. Use `--study-preset` only with the study's material assignments.

Figures 1 and 2 have four equal panel areas, external bold panel letters and no panel titles. Their image aspect ratios are preserved. Figure 2 discovers the available spectrum keys and marks non-selection without an REA box or selection line. A selected-size square is schematic, not the location of the actual spectral window. Spectra are peak-normalized for display only; selection uses unnormalized spectra. `--pixel-size-mm` converts support and wavenumber axes assuming `dx=1`; image coordinates remain in pixels.

Figure 3 uses all supplied map sets and response-image IDs. Without a phase table it contains three property-response panels; with a phase table it also includes composition. Reduced image populations are identified from the actual response data. Single-window SEM remains undefined. The assembled outputs are `Figure_1.png/.pdf`, `Figure_2.png/.pdf` and `Figure_3.png/.pdf`, with no panel titles and bold panel labels. All three scripts default to `outputs/figures/`; Figure 3 also writes numerical summaries under `derived_data/` within its output directory. Use `--dpi` for lower-resolution previews.

## Read and interpret the results

These reference values describe the dune-sand dataset, not expected results for new images.

### Structural support and numerical conventions

The study's baseline protocol gives `L_REA = 1536` pixels, or approximately **2.01 mm**, against a 2048-pixel finite reference. The screening masks retain approximately 15.23%, 39.75%, 37.35%, 61.48%, 39.39%, 45.35% and 40.69% of BSE images 1–7, respectively.

The mask selects admissible **whole windows**. It is not multiplied into the image, and windows may contain up to 2% outside-mask pixels under the stated coverage threshold. Material outside the selected domain remains part of the specimen; its response is not established by an analysis restricted to admissible windows. Local marginal compatibility also does not establish full spatial stationarity.

Several settings have distinct roles:

- The requested **121 sampling positions** form an 11 by 11 grid, with the centre included if distinct. A separate preliminary availability check in the core program probes **49 positions**.
- The mean residual used for selection undergoes five-point median regularization and a suffix-maximum envelope. The optional three-point plotted guides serve only display purposes.
- The structural selection combines the regularized mean residual and the unsmoothed spectral residual. Ensemble variability is an auxiliary diagnostic and does not enter the acceptance rule.
- The persistent-tail rule requires at least two nonreference supports. With the supplied support sequence, 1536 pixels is the largest selectable side; the upper end of a selection distribution is therefore constrained by the tested grid.

### Uncertainty and sensitivity results

| Analysis | Result for the supplied dataset |
| --- | --- |
| Baseline protocol | 1536 pixels. |
| Whole-image bootstrap, 10,000 draws, seed 20261029 | 3310 draws yield no selected support: 33.10%. |
| Successful bootstrap selections only | 2.5–97.5 percentile range: 896–1536 pixels. |
| Masked sampling with 49 requested positions | No selected support. |
| Unmasked sampling with 121 requested positions | 1280 pixels, for a different sampling population. |
| Screening-window side 1536 pixels | Only two common valid supports remain; no estimate. |
| Screening-window side 3072 pixels | 1536 pixels. |
| Reference side 1536 / 1792 pixels | No selection / 1280 pixels. |
| Reference sides 2560–4096 pixels | At least one image fails the production availability check. |

The successful-only bootstrap range is **not an unconditional 95% confidence interval**. The bootstrap resamples whole images while retaining their dependence across support sizes; masks and numerical settings are fixed. It assumes that the seven image fields can be treated as independent realizations and does not include uncertainty in resolution, segmentation, calibration or mask thresholds.

The reference support is only 1.33 times the selected side and has not been demonstrated to reach an asymptotic limit. The results therefore characterize a conditional support under the stated finite-reference protocol. Pixel size and resolved microstructure may affect the estimate. Validation on this material and two-dimensional geometry does not establish a universal support for other materials or a bulk three-dimensional REV.

The supplementary reference-sensitivity file also records a **2560-pixel exploratory calculation** that uses 121-position eligibility instead of the production 49-position gate. Its `enforces_production49_common_size_gate` flag is false. Keep that calculation separate from results produced under the production eligibility rule.

### Property responses

The property-map screening leaves admissible windows in images **1, 5 and 6**. The listed support sizes from 32 through 320 pixels have 12 pooled windows each, with four per image. The 384-pixel result contains only four windows from image 6; there are no accepted windows at 448 or 512 pixels. The dataset contains **160 solved windows** in total.

The axial quantities `Ex = 1/S11` and `Ey = 1/S22` are derived from the full plane-strain compliance. The shear quantity `Gxy = 1/S33` retains coupling through that matrix inversion and is generally different from `C33`. The pooled-window standard errors describe the available windows; they are not confidence intervals based on independent images. The image-6-only result should be interpreted with its different sampling population.

Recorded cell-solve times total approximately 420.79 seconds. Hardware and thread settings for these recorded timings are not recorded, so this value is not a runtime guarantee for another machine.

### File formats

CSV files contain numerical tables. NPZ files contain NumPy arrays and can be opened with `numpy.load(path, allow_pickle=False)`. Checkpoints use gzip-compressed JSON with explicit array shapes and dtypes. In Figure 2 NPZ files, numeric NaN in `accepted_L` marks non-selection. In result JSON, `null` identifies an undefined or unavailable value according to the field. In `bootstrap_draws.npz`, an infinite value in `selected_pixels` is a non-selection marker, not a measurement of an infinite physical support.

## Troubleshooting

| Symptom | Action |
| --- | --- |
| No images or missing masks | Check the glob/extension, unique image stems, mask suffix and matching dimensions. |
| Checkpoint dataset/configuration mismatch | Start a new checkpoint directory for changed inputs or settings; do not mix records from different runs. |
| Incomplete checkpoint set | Process every input position before aggregation; retain the same complete list for individual jobs. |
| Unexpected seven-image result | Pass your own `--checkpoint-dir` to every aggregation/reconstruction command and your new `--data-file` to Figure 2. |
| No common support or no selected REA | Inspect mask coverage, image dimensions, support grid, reference and tolerances. Report non-selection and sensitivity; a finite REA is not guaranteed. |
| No property windows at large supports | Inspect per-image availability. The response population may shrink with increasing support. |
| Memory exhaustion | Reduce input dimensions only with a justified resolution study, use a larger-memory machine, or use the sequential checkpoint workflow for structural uncertainty. |
| Wrong phase labels on a new material | Omit `--study-preset` and supply your own phase table only when its assignments are valid. |

## Reproducibility, citation and license

The inputs underlying the reference statistics are identified by [Zenodo version DOI 10.5281/zenodo.21717613](https://doi.org/10.5281/zenodo.21717613). The [concept DOI](https://doi.org/10.5281/zenodo.21717612) identifies the dataset's version series.

When using the repository, cite the accompanying study, the **specific Zenodo version used**, and the **Git commit or release used**. Record the command-line settings and any dataset-specific changes. The bundled numerical check compares scientific values with tolerances and excludes elapsed times and source-hash metadata from numerical equality tests; it verifies core source files separately. Checkpoint formats and source verification are documented in [validation/README.md](validation/README.md).

The software is distributed under the GNU General Public License v3; see [LICENSE](LICENSE). Data licensing is specified in the selected Zenodo record.

For questions about the method or data, contact **Dr. Christian Tantardini** at [christiantantardini@ymail.com](mailto:christiantantardini@ymail.com) or **Prof. Eduardo Garzanti** at [eduardo.garzanti@unimib.it](mailto:eduardo.garzanti@unimib.it). Reproducible software problems can be reported through [GitHub Issues](https://github.com/Christian48596/representativeness-emerges/issues), including the command, software environment and error message.
