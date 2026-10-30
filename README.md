# Statistical admissibility and long-wavelength structural convergence determine representative support in granular materials

**Christian Tantardini, Fernando Alonso-Marroquín, Abdullah Alqubalee and Eduardo Garzanti**

This repository provides a Python workflow for determining the spatial support over which the structure of a granular material can be represented consistently and for testing that support against independent physical properties. It combines local statistical screening, convergence of image means and long-wavelength spatial correlations, and periodic thermal and elastic homogenization.

The accompanying study applies this framework to two-dimensional polished sections of dune sand. Backscattered-electron (BSE) images provide the structural information. Numerical property maps derived from QEMSCAN mineralogical classifications provide independent conductivity and elasticity fields. The repository contains the analysis programs, reference statistics, uncertainty and sensitivity calculations, and scripts for producing the manuscript figures.

## Contents

- [Scientific workflow](#scientific-workflow)
- [Repository organization](#repository-organization)
- [Installation and a first verification](#installation-and-a-first-verification)
- [Download and prepare the data](#download-and-prepare-the-data)
- [Run the complete analysis](#run-the-complete-analysis)
- [Quantify uncertainty and parameter sensitivity](#quantify-uncertainty-and-parameter-sensitivity)
- [Generate the figures](#generate-the-figures)
- [Read and interpret the results](#read-and-interpret-the-results)
- [Apply the programs to another dataset](#apply-the-programs-to-another-dataset)
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
| [`rea_stationary_mask_cli_v24_general.py`](rea_stationary_mask_cli_v24_general.py) | Evaluate structural statistics and spectra across support sizes and select a persistent structural support. |
| [`property_full_homogenized_inside_property_masks.py`](property_full_homogenized_inside_property_masks.py) | Solve periodic thermal and elastic cell problems inside admissible property-map domains. |
| [`validation/`](validation/) | Bootstrap, sensitivity, numerical reconstruction and verification scripts. |
| [`validation/checkpoints/`](validation/checkpoints/) | Compressed per-image statistics and spectra for reproducing the uncertainty analysis without processing the large images again. |
| [`validation/results/`](validation/results/) | Reference bootstrap, sensitivity and per-image results. |
| [`validation/inputs/`](validation/inputs/) | Archived structural summary used for numerical comparison. |
| [`figures/`](figures/) | Scripts for Figures 1–3 and the supplementary property-map illustration. |
| [`derived_data/`](derived_data/) | Compact numerical plotting inputs, phase fractions and property summaries. |
| [`docs/ORIGINAL_RUN_COMMANDS.md`](docs/ORIGINAL_RUN_COMMANDS.md) | Additional documentation of the study's computational settings. |
| [`requirements.txt`](requirements.txt) | Required Python packages. |
| [`requirements-revision.txt`](requirements-revision.txt) | Pinned package versions used to verify the numerical reconstruction. |

Large images, full-resolution masks and numerical property maps are distributed through Zenodo. Generated outputs are written under `outputs/` by the commands below.

## Installation and a first verification

The workflow uses NumPy, SciPy, Matplotlib and Pillow. The recorded verification environment uses Python 3.12.14. The shell examples below are written for **Bash on Linux or macOS**; use a Bash terminal when copying blocks containing arrays.

Clone the repository and create a virtual environment:

```bash
git clone https://github.com/Christian48596/representativeness-emerges.git
cd representativeness-emerges

python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

To use the recorded package versions instead, install them in a Python 3.12 environment:

```bash
python -m pip install -r requirements-revision.txt
```

A first numerical verification requires only the files in this repository:

```bash
python validation/check_reference_results.py
```

This command reconstructs the baseline, all 10,000 bootstrap draws, the leave-one-image-out results, the sampling and screening sensitivity analyses, and the numerical inputs for Figure 2. It compares them with the bundled reference results and checks the three core program files against their recorded Git blob identifiers. A successful run prints `PASS`, followed by the baseline support and bootstrap summary.

This compact check uses previously calculated per-image statistics. It verifies the numerical reconstruction without repeating full-resolution image processing or finite-element solves. The next sections explain how to perform those calculations from the raw data.

For the command-line interfaces:

```bash
python detect_stationary_domains_general.py --help
python rea_stationary_mask_cli_v24_general.py --help
python property_full_homogenized_inside_property_masks.py --help
```

Full-resolution processing is memory intensive: each BSE field contains approximately 205 million pixels, and intermediate arrays and sampled windows require additional memory. The compact verification is the appropriate starting point on machines with limited memory.

## Download and prepare the data

The dataset is available through [Zenodo](https://doi.org/10.5281/zenodo.21717612). Download the following files from one selected version of the record:

| File | Contents |
| --- | --- |
| `BSE_full_resolution_images.tar.gz` | Seven full-resolution BSE images. |
| `QEMSCAN_numerical_property_maps.tar.gz` | Numerical conductivity, Young's-modulus and Poisson-ratio maps. |
| `reference_outputs_and_masks.tar.gz` | Admissibility masks, numerical reference outputs and execution logs; additional contents depend on the selected record version. |
| `README_DATA.md` | Dataset documentation. |
| `SHA256SUMS.txt` | Checksums for that record version. |

Create a `data/` directory in the repository and place the five downloaded files there. From that directory, verify the checksums with `sha256sum -c SHA256SUMS.txt` on Linux or `shasum -a 256 -c SHA256SUMS.txt` on macOS.

From the repository root, extract the archives:

```bash
mkdir -p data

tar -xzf data/BSE_full_resolution_images.tar.gz -C data
tar -xzf data/QEMSCAN_numerical_property_maps.tar.gz -C data
tar -xzf data/reference_outputs_and_masks.tar.gz -C data

export REA_DATA_ROOT="$PWD/data"
```

The extracted data directory must contain:

| Relative path | Meaning |
| --- | --- |
| `BSE/1_cut.tif` … `BSE/7_cut.tif` | BSE image fields. |
| `QEMSCAN/kappaMap_1.csv` … `kappaMap_7.csv` | Thermal conductivity, in W m^-1 K^-1. |
| `QEMSCAN/EMap_1.csv` … `EMap_7.csv` | Young's modulus, in GPa. |
| `QEMSCAN/nuMap_1.csv` … `nuMap_7.csv` | Poisson ratio, dimensionless. |
| `reference_outputs/` | Archived masks, summaries and logs. |

The BSE pixel size used for physical length reporting is **1.31 micrometres**. The QEMSCAN property-map pixel size is **0.01 mm**. The structural spectrum is calculated in pixel coordinates using `--pixel-size 1`; physical side lengths are reported through a separate conversion.

The numerical QEMSCAN CSVs are already classified and assigned material properties. The repository starts from these arrays; it does not reconstruct the upstream mineral classification from RGB illustrations. Use the CSV values directly. Intensity-rescaled visualization images do not preserve the physical property values required by the solver. The three property maps for each image must have matching dimensions and registration.

## Run the complete analysis

Run the following blocks in sequence from the repository root, with the virtual environment active. Check that each stage finishes successfully before starting the next one. Commands run in the foreground and write detailed progress to log files.

### Prepare the input lists and a separate calculation directory

The following setup writes fresh results to `outputs/full_run/reference_outputs/`. Symbolic links give the run directory the same layout as the downloaded dataset, allowing the plotting and validation scripts to use the fresh results later.

```bash
REA_INPUT_ROOT="$REA_DATA_ROOT"
REA_RUN_ROOT="$PWD/outputs/full_run"
REA_RUN_DIR="$REA_RUN_ROOT/reference_outputs"
mkdir -p "$REA_RUN_DIR"

if [ ! -e "$REA_RUN_ROOT/BSE" ]; then
  ln -s "$REA_INPUT_ROOT/BSE" "$REA_RUN_ROOT/BSE"
fi
if [ ! -e "$REA_RUN_ROOT/QEMSCAN" ]; then
  ln -s "$REA_INPUT_ROOT/QEMSCAN" "$REA_RUN_ROOT/QEMSCAN"
fi

bse_files=()
kappa_files=()
elastic_files=()
poisson_files=()
for i in {1..7}; do
  bse_files+=("$REA_INPUT_ROOT/BSE/${i}_cut.tif")
  kappa_files+=("$REA_INPUT_ROOT/QEMSCAN/kappaMap_${i}.csv")
  elastic_files+=("$REA_INPUT_ROOT/QEMSCAN/EMap_${i}.csv")
  poisson_files+=("$REA_INPUT_ROOT/QEMSCAN/nuMap_${i}.csv")
done
```

Keep this Bash session open so the variables and arrays remain available. The three core programs receive their inputs explicitly; the additional validation and figure scripts also recognize `REA_DATA_ROOT`.

### 1. Detect admissible domains in the BSE images

```bash
python detect_stationary_domains_general.py "${bse_files[@]}" \
  --input-type image --field-transform raw --max-side 0 \
  --output-dir "$REA_RUN_DIR/stationary_domains_FULLRES_R2048_STRICT" \
  --stationarity-window 2048 --stride 512 \
  --tau-mu 0.02 --tau-sigma 0.06 --min-component-windows 4 --dpi 300 \
  > "$REA_RUN_DIR/run_detect_stationary_domains_FULLRES_R2048_STRICT.log" 2>&1
```

This stage evaluates local means and standard deviations on moving 2048-pixel windows, with a stride of 512 pixels. It applies the configured tolerances, morphological opening/closing and connected-component filtering, then exports the pixel masks and local diagnostics.

The mask files are named `1_cut_stationary_mask.tif` through `7_cut_stationary_mask.tif`. Their values indicate the domain available for subsequent window selection. The screening tests local marginal compatibility; the interpretation of that domain is discussed under [Read and interpret the results](#read-and-interpret-the-results).

### 2. Determine the structural support

```bash
python rea_stationary_mask_cli_v24_general.py "${bse_files[@]}" \
  --output-dir "$REA_RUN_DIR/combined_REA_final_STRICT_v24" \
  --output-prefix REA_stationary_domains_STRICT_v24 \
  --representative-image 1 --max-side 0 --pixel-size 1 --field-mode raw \
  --stationary-mask-dir "$REA_RUN_DIR/stationary_domains_FULLRES_R2048_STRICT" \
  --stationary-mask-suffix _stationary_mask.tif --min-mask-fraction 0.98 \
  --window-sizes 256,384,512,640,768,896,1024,1280,1536,1792,2048 \
  --reference-L 2048 --selected-spectra 256,512,1024,2048 \
  --n-common-k-grid 600 --n-window-samples 121 --errorbar-mode sem \
  --spectral-band low --low-k-fraction 0.25 --low-k-weight inverse \
  --tau-Z 0.025 --tau-C 0.30 --tau-ens 0.10 --tail-fraction 0.5 \
  --panel-b-mode raw+smooth --panel-b-smooth-window 3 \
  --panel-c-mode raw+smooth --panel-c-smooth-window 3 --dpi 600 --no-pdf \
  > "$REA_RUN_DIR/run_REA_stationary_STRICT_v24.log" 2>&1
```

For each support, candidate windows must have at least 98% mask coverage. The sampled means characterize their gray-level average, while a selected admissible window from each image supplies the spatial spectrum. The ensemble analysis compares support-dependent statistics with the 2048-pixel reference and applies the persistent joint selection rule.

For the supplied dataset, the baseline selection is **1536 BSE pixels**, corresponding to **2.01216 mm**. The principal numerical output is:

```text
outputs/full_run/reference_outputs/combined_REA_final_STRICT_v24/REA_stationary_domains_STRICT_v24_combined_window_statistics.csv
```

Use a basename for `--output-prefix`: the program joins it to `--output-dir`. The explicit settings above reproduce the study's full-resolution protocol; running the program with its defaults selects different preprocessing and sampling settings.

### 3. Detect admissible domains in the property maps

```bash
python detect_stationary_domains_general.py \
  --kappa-maps "${kappa_files[@]}" \
  --E-maps "${elastic_files[@]}" \
  --nu-maps "${poisson_files[@]}" \
  --stationarity-source combined --combined-transform log-robust --max-side 0 \
  --output-dir "$REA_RUN_DIR/stationary_property_combined_STRICT" \
  --stationarity-window 256 --stride 64 \
  --tau-mu 0.02 --tau-sigma 0.06 --min-component-windows 4 --dpi 300 \
  > "$REA_RUN_DIR/run_detect_stationary_property_combined_STRICT.log" 2>&1
```

This stage forms a combined screening field from the numerical conductivity, Young's-modulus and Poisson-ratio maps. The conductivity and modulus channels are log transformed and robustly standardized; Poisson ratio is robustly standardized without a logarithmic transform. Masks are generated on the property-map grid and are named `property_combined_1_stationary_mask.tif` through `property_combined_7_stationary_mask.tif`.

### 4. Solve the periodic property problems

```bash
python property_full_homogenized_inside_property_masks.py \
  --kappa-maps "${kappa_files[@]}" \
  --E-maps "${elastic_files[@]}" \
  --nu-maps "${poisson_files[@]}" \
  --mask-dir "$REA_RUN_DIR/stationary_property_combined_STRICT" \
  --mask-template 'property_combined_{i}_stationary_mask.tif' \
  --output-dir "$REA_RUN_DIR/property_full_homogenized_property_masks_STRICT" \
  --output-prefix figure5_property_full_homogenized_property_masks_STRICT \
  --window-sizes 32,48,64,80,96,112,128,160,192,204,224,256,320,384,448,512 \
  --L-rea 204 --property-pixel-size-mm 0.01 --min-mask-fraction 0.98 \
  --n-window-samples 4 --elastic-mode plane_strain --errorbar-mode sem \
  --rng-seed 12345 --dpi 600 --no-pdf \
  > "$REA_RUN_DIR/run_property_full_homogenized_property_masks_STRICT.log" 2>&1
```

The nominal transferred support is **204 property-map pixels**, or **2.04 mm**. The solver evaluates two thermal and three elastic loading cases using periodic Q4 finite elements, 2 by 2 Gauss integration, zero-mean constraints and a sparse direct solve.

The output directory contains a per-window response CSV, a support-wise summary CSV and response curves. The per-window file is `figure5_property_full_homogenized_property_masks_STRICT_window_values.csv`; it records conductivity, stiffness, axial moduli and elapsed times. Figure postprocessing inverts the full symmetric 3 by 3 stiffness matrix to obtain the compliance and shear modulus, retaining normal–shear coupling.

After all four stages have succeeded, point the auxiliary scripts to this complete fresh run:

```bash
export REA_DATA_ROOT="$REA_RUN_ROOT"
```

To work with the downloaded reference outputs instead, use `export REA_DATA_ROOT="$PWD/data"`. This choice determines which images, masks and cell responses the subsequent data-dependent commands read.

## Quantify uncertainty and parameter sensitivity

The validation scripts assess how image sampling, candidate positions, mask screening and finite-reference choice affect the selected support. They use the core spectral, ensemble-aggregation and selection routines.

### Reconstruct the stored analyses without raw-image processing

From the repository root:

```bash
python validation/rea_validation_wrapper.py --aggregate --bootstrap 10000
python validation/rstat_sensitivity_wrapper.py --aggregate
python validation/rebuild_derived.py
```

These commands read `validation/checkpoints/` and write new numerical files to `outputs/validation/`. They reconstruct the stored analyses even when `REA_DATA_ROOT` points to a fresh run: **the default checkpoint source remains the bundled directory**.

| Generated file | Contents |
| --- | --- |
| `validation_results.json` | Baseline, masked and unmasked sampling comparisons, reference-availability checks, whole-image bootstrap and leave-one-image-out results. |
| `rstat_sensitivity_results.json` | Results for screening-window sides 1536 and 3072 pixels at fixed stride. |
| `reference_sensitivity_extra.json` | Shorter-reference tests and explicitly flagged exploratory reference calculations. |
| `bootstrap_draws.npz` | Image indices, selected supports and spectral results for every bootstrap draw. |
| `baseline_per_image_1.csv` … `baseline_per_image_7.csv` | Individual-image statistics across support sizes. |
| `figure2_plot_data.npz` | Numerical series and spectra used by the Figure 2 plotter. |

### Recompute the image statistics and screening masks

To recompute the underlying statistics from the images and the selected data root, run the image jobs sequentially and explicitly pass their checkpoint directory to every aggregation command:

```bash
for i in 1 2 3 4 5 6 7; do
  python validation/rea_validation_wrapper.py \
    --image "$i" --data-root "$REA_DATA_ROOT" \
    --checkpoint-dir outputs/validation/checkpoints
  python validation/rstat_sensitivity_wrapper.py \
    --image "$i" --data-root "$REA_DATA_ROOT" \
    --checkpoint-dir outputs/validation/checkpoints
done

python validation/rea_validation_wrapper.py \
  --aggregate --bootstrap 10000 --data-root "$REA_DATA_ROOT" \
  --checkpoint-dir outputs/validation/checkpoints
python validation/rstat_sensitivity_wrapper.py \
  --aggregate --data-root "$REA_DATA_ROOT" \
  --checkpoint-dir outputs/validation/checkpoints
python validation/rebuild_derived.py \
  --checkpoint-dir outputs/validation/checkpoints
```

The structural wrapper reads the baseline masks under `REA_DATA_ROOT/reference_outputs/stationary_domains_FULLRES_R2048_STRICT/`. The screening-sensitivity wrapper recomputes alternative masks at 1536 and 3072 pixels. `--output-dir` changes the numerical-output destination; `--checkpoint-dir` controls which per-image records are written or read. `rebuild_derived.py` reads checkpoints and does not take a `--data-root` argument.

The reconstruction wrappers encode the supplied study's seven-image protocol and baseline checks. They are not a general bootstrap interface for an arbitrary number of images. See [validation/README.md](validation/README.md) for the checkpoint format and verification details.

## Generate the figures

The figure scripts use numerical arrays and archived or freshly computed outputs. They write PNG and PDF files.

| Figure | Scientific content | Main inputs |
| --- | --- | --- |
| Figure 1 | Image field, admissible domain and local screening diagnostics | BSE image 4, its mask and window-diagnostic CSV. |
| Figure 2 | Structural mean, spectral convergence and selected spectra | BSE image 1 and a numerical plotting NPZ. |
| Figure 3 | Phase fractions, conductivity, stiffness and axial/shear moduli | Seven property maps and solved cell responses, or the equivalent compact inputs. |
| Supplementary property-map figure | Mineral/property-map illustration | Numerical maps for image 1. |

### Figures from the downloaded reference data

```bash
python figures/regenerate_figure1_layout.py --data-root "$PWD/data"
python figures/regenerate_figure2_layout.py --data-root "$PWD/data"
python figures/regenerate_revised_figures.py --data-root "$PWD/data"
```

Figure 2 uses the bundled `derived_data/figure2_plot_data.npz` by default. Figures 1 and 2 are written under `outputs/figures/`; Figure 3 and the supplementary map figure are written under `outputs/revised_figures/IMG/`, with derived tables under `outputs/revised_figures/derived_data/`.

Figure 3 and the supplementary map illustration can also be generated without downloading the raw data, using only the compact inputs:

```bash
env -u REA_DATA_ROOT python figures/regenerate_revised_figures.py
```

Unsetting the environment variable for this one command ensures that the plotter selects its bundled inputs.

### Figures from a fresh full calculation

After completing the four-stage calculation and recomputing the validation checkpoints above:

```bash
python figures/regenerate_figure1_layout.py --data-root "$REA_RUN_ROOT"
python figures/regenerate_figure2_layout.py \
  --data-root "$REA_RUN_ROOT" \
  --data-file outputs/validation/figure2_plot_data.npz
python figures/regenerate_revised_figures.py --data-root "$REA_RUN_ROOT"
```

The explicit `--data-file` directs Figure 2 to the newly reconstructed curves. Each plotter accepts `--output-dir` if separate destinations are needed for archived and fresh figures. The programs retain their existing output filenames, including the `_revised` suffix.

Figures 1 and 2 use four equal plot areas and external bold panel letters. Spectral curves in Figure 2 are peak-normalized for display; the selection uses the unnormalized spectra. Its white square indicates the selected support size and does not identify the actual admissible spectral-window location. Figure 3 uses a linear phase-fraction scale and derives `Gxy = 1/S33` from the full compliance matrix. Further plotting details are in [figures/README.md](figures/README.md).

## Read and interpret the results

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

The supplementary reference-sensitivity file also records a **2560-pixel exploratory calculation** that uses 121-position eligibility instead of the production 49-position gate. Its `enforces_original49_common_size_gate` flag is false. Keep that calculation separate from results produced under the production eligibility rule.

### Property responses

The property-map screening leaves admissible windows in images **1, 5 and 6**. The listed support sizes from 32 through 320 pixels have 12 pooled windows each, with four per image. The 384-pixel result contains only four windows from image 6; there are no accepted windows at 448 or 512 pixels. The dataset contains **160 solved windows** in total.

The axial quantities `Ex = 1/S11` and `Ey = 1/S22` are derived from the full plane-strain compliance. The shear quantity `Gxy = 1/S33` retains coupling through that matrix inversion and is generally different from `C33`. The pooled-window standard errors describe the available windows; they are not confidence intervals based on independent images. The image-6-only result should be interpreted with its different sampling population.

Recorded cell-solve times total approximately 420.79 seconds. Hardware and thread settings for those historical timings are not recorded, so this value is not a runtime guarantee for another machine.

### File formats

CSV files contain numerical tables. NPZ files contain NumPy arrays and can be opened with `numpy.load(path, allow_pickle=False)`. Checkpoints use gzip-compressed JSON with explicit array shapes and dtypes. In result JSON, `null` identifies an undefined or unavailable value according to the field. In `bootstrap_draws.npz`, an infinite value in `selected_pixels` is a non-selection marker, not a measurement of an infinite physical support.

## Apply the programs to another dataset

Use the three core programs' command-line options to supply new fields and choose analysis parameters. Define the pixel calibration, candidate support sequence, reference size, screening window and tolerances for the material and imaging conditions being studied. Check that admissible supports exist across the image ensemble before interpreting a selected size.

For property calculations, supply registered numerical conductivity, Young's-modulus and Poisson-ratio maps in consistent units, with masks defined on the same grid. Choose the appropriate mechanical assumption and document its relationship to the physical problem. A structural support determined from another image channel should be assessed against independent property responses at the corresponding physical length.

The figures, phase assignments and seven-image validation wrappers in this repository are configured for the dune-sand study. Adapt those dataset-specific components before using them elsewhere. Report non-selections and sensitivity to the reference, screening and sampling choices along with any selected support.

## Reproducibility, citation and license

The three core analysis programs correspond to GitHub commit [`8bd49ef71bbd8e320bdb00b3633f8d137cb72bb0`](https://github.com/Christian48596/representativeness-emerges/tree/8bd49ef71bbd8e320bdb00b3633f8d137cb72bb0). The archived inputs underlying the reference statistics are identified by [Zenodo version DOI 10.5281/zenodo.21717613](https://doi.org/10.5281/zenodo.21717613). The [concept DOI](https://doi.org/10.5281/zenodo.21717612) identifies the dataset's version series.

When using the repository, cite the accompanying study, the **specific Zenodo version used**, and the **Git commit or release used**. Record the command-line settings and any dataset-specific changes. The bundled numerical check compares scientific values with tolerances and excludes elapsed times and source-hash metadata from numerical equality tests; it verifies core source files separately. Historical source-hash details and checkpoint provenance are documented in [validation/README.md](validation/README.md).

The software is distributed under the GNU General Public License v3; see [LICENSE](LICENSE). Data licensing is specified in the selected Zenodo record.

For questions about the method or data, contact **Dr. Christian Tantardini** at [christiantantardini@ymail.com](mailto:christiantantardini@ymail.com) or **Prof. Eduardo Garzanti** at [eduardo.garzanti@unimib.it](mailto:eduardo.garzanti@unimib.it). Reproducible software problems can be reported through [GitHub Issues](https://github.com/Christian48596/representativeness-emerges/issues), including the command, software environment and error message.
