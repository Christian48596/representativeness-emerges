# Structural uncertainty analysis for N images

These scripts calculate whole-image bootstrap uncertainty, leave-one-image-out selections, sampling-grid comparisons, screening-window sensitivity and reference-support diagnostics. They import the repository's original numerical routines. Image count and filenames are supplied by the dataset; there is no seven-image limit.

## Verify the published study

From the repository root:

```bash
python validation/check_reference_results.py
```

This is specifically the regression check for the archived seven-image study. It uses bundled checkpoints, repeats all 10,000 bootstrap draws, compares every recorded scientific value and Figure 2 array, and checks the three core source files against their Git blob identifiers. It does not analyse a user's new data. Recorded wall times and historical source-hash metadata are excluded from numerical comparison.

## Process another dataset

Use an ordered image list and either a matching ordered mask list (`--masks`) or a mask directory. Directory-based mask names are the image stem followed by `_stationary_mask.tif`; `--mask-suffix` changes the suffix. Images may have arbitrary names. The mask and its image must have identical dimensions. Masks should be produced by the domain detector before running the structural wrapper.

The following Bash example processes all supplied images sequentially and stores independent outputs. Adjust the support and screening scales to the image dimensions and the physical problem. Use a fresh output/checkpoint directory for a different dataset or configuration.

```bash
COMMON_ARGS=(
  --images /absolute/path/to/images/*.tif
  --mask-dir /absolute/path/to/masks
  --window-sizes 64,96,128,192,256
  --reference-L 256
  --extra-reference-sizes 320,384
  --n-window-samples 121
  --n-common-k-grid 600
  --tau-Z 0.025 --tau-C 0.30
  --min-mask-fraction 0.98 --tail-fraction 0.5
  --screening-sizes 128,256 --screening-stride 32
  --screening-tau-mu 0.02 --screening-tau-sigma 0.06
  --min-component-windows 4
  --checkpoint-dir outputs/my_data/checkpoints
  --output-dir outputs/my_data
)

python validation/rea_validation_wrapper.py --all-images "${COMMON_ARGS[@]}"
python validation/rstat_sensitivity_wrapper.py --all-images "${COMMON_ARGS[@]}"

python validation/rea_validation_wrapper.py --aggregate \
  --checkpoint-dir outputs/my_data/checkpoints \
  --output-dir outputs/my_data --bootstrap 10000 --seed 20261029
python validation/rstat_sensitivity_wrapper.py --aggregate \
  --checkpoint-dir outputs/my_data/checkpoints --output-dir outputs/my_data
python validation/rebuild_derived.py \
  --checkpoint-dir outputs/my_data/checkpoints --output-dir outputs/my_data
```

`--all-images` uses every supplied path, so the same commands work for nine, twenty or another positive number of images. `--image INDEX` instead processes one 1-based position in that complete list. Always supply the complete image list when producing individual checkpoints; this lets aggregation identify incomplete runs.

Without `--images`, raw processing discovers supported image files under `REA_DATA_ROOT/BSE/` in natural filename order. `--data-root` overrides the environment variable. The default baseline mask directory is `reference_outputs/stationary_domains_FULLRES_R2048_STRICT/` below the data root. Explicit image and mask paths are independent of that directory convention.

## Checkpoint configuration and consistency

Each newly calculated checkpoint includes the complete ordered image/mask list, the image and mask file sizes and modification timestamps, and the numerical configuration. A `dataset_manifest.json` file identifies the expected ensemble. Aggregation infers N and settings from this metadata; it does not require the original images to be available when reading a complete checkpoint set.

Aggregation rejects missing IDs, incomplete image sets, duplicate IDs, mixed configurations and mismatched manifests. Raw writes reject an existing checkpoint directory when the image list, file identities or configuration differ. Use a new directory when changing input files or analysis settings. `--expected-images N` provides an additional explicit image-count check.

Use the same `COMMON_ARGS` for structural and screening-window image jobs, because their shared manifest records the complete configuration. On aggregation, settings are inferred automatically. Explicit settings that contradict stored values are rejected. This prevents results calculated at different scales or thresholds from being silently pooled.

For new image jobs the default checkpoint destination is `outputs/validation/checkpoints/`. For aggregation without `--checkpoint-dir`, the default input is the bundled `validation/checkpoints/`. Therefore, always provide `--checkpoint-dir` for a custom dataset.

The compact bundled study can be reconstructed without raw images using:

```bash
python validation/rea_validation_wrapper.py --aggregate --bootstrap 10000
python validation/rstat_sensitivity_wrapper.py --aggregate
python validation/rebuild_derived.py
```

The general aggregation command does not assert a 1536-pixel result or compare arbitrary data with the dune-sand reference. `--verify-reference` explicitly enables the seven-image study check; `check_reference_results.py` uses that option and the bundled original CSV.

## Parameters

| Option | Meaning |
| --- | --- |
| `--window-sizes` | Comma-separated support sides in pixels; at least four distinct sizes. |
| `--reference-L` | Requested reference side, included in the support sequence. |
| `--extra-reference-sizes` | Additional support sides for reference sensitivity; an empty string disables them. |
| `--n-window-samples` | Requested candidate sampling positions; the production availability gate remains 49 positions. |
| `--n-common-k-grid` | Number of points in the common radial wavenumber grid. |
| `--tau-Z`, `--tau-C` | Positive mean and spectral tolerances. |
| `--min-mask-fraction` | Required coverage of an analysis window by the mask. |
| `--tail-fraction` | Fraction used by the persistent-tail selection rule. |
| `--low-k-fraction`, `--low-k-weight` | Low-wavenumber band and `uniform` or `inverse` weighting. |
| `--screening-sizes`, `--screening-stride` | Alternative screening-window sides and stride, in pixels. |
| `--screening-tau-mu`, `--screening-tau-sigma` | Local screening tolerances for the alternative masks. |
| `--min-component-windows` | Minimum connected coarse-grid component size. |
| `--short-reference-sizes` | Shorter references for reconstruction; unavailable values are omitted or reported as unavailable. |
| `--selected-spectra` | Support sides included in the Figure 2 plotting archive. |
| `--pixel-size-mm` | Optional physical reporting calibration. Spectral calculations retain `dx=1` in pixel coordinates. |
| `--bootstrap`, `--seed` | Bootstrap count and random seed for structural aggregation. |

Without explicit settings, the recorded study defaults apply: supports 256,384,512,640,768,896,1024,1280,1536,1792,2048; reference 2048; 121 requested positions; 600 radial points; mean/spectral tolerances 0.025/0.30; minimum mask coverage 0.98; tail fraction 0.5; low-wavenumber fraction 0.25 with inverse weighting. The screening sensitivity defaults to 1536/3072-pixel windows, stride 512, tolerances 0.02/0.06 and four connected coarse cells. These defaults are not recommended universal scales for another dataset.

## Outputs and interpretation

`validation_results.json` contains the baseline, grid/unmasked comparisons, larger-reference eligibility, bootstrap and leave-one-image-out results. `bootstrap_draws.npz` stores a `(number_of_draws, N)` index array, selected supports and spectral results. `rstat_sensitivity_results.json` records the alternative screening calculations. `reference_sensitivity_extra.json` records shorter-reference checks and separately flagged exploratory larger-reference calculations.

The common-support gate is applied across the complete ensemble. If fewer than four supports remain, the baseline reports no estimate with a reason. If the requested reference is unavailable but enough common supports exist, the largest valid support is used, as in the core program; the baseline records the requested reference when it differs from the effective one. No selection is a valid reported outcome.

Bootstrap draws resample whole-image records, preserving their within-image dependence across support sizes. They use fixed masks, settings and the common support grid. The result assumes that the supplied image fields represent independent realizations; changing N alone does not establish that assumption. Successful-only percentiles remain conditional when some draws have no selection.

With N=1, image bootstrap resampling is degenerate. The result explicitly marks that independent-image uncertainty cannot be estimated and leaves the leave-one-image-out list empty. With no common admissible supports, outputs retain unavailable values and the plotting archive uses NaN curves and NaN selected support.

Result JSON uses `null` for undefined values. In `bootstrap_draws.npz`, `selected_pixels=inf` marks non-selection, not an infinite measured support. Figure 2 NPZ arrays use numeric `NaN` for no selection and support `numpy.load(..., allow_pickle=False)`. Pass the newly generated plotting archive explicitly:

```bash
python figures/regenerate_figure2_layout.py \
  --image-file /absolute/path/to/a/representative_image.tif \
  --data-file outputs/my_data/figure2_plot_data.npz
```

The larger-reference exploratory branch uses the actual sampling grid without enforcing the production 49-position gate and is labelled accordingly. It must not be interpreted as a result of the production eligibility rule.

## Formats and provenance

Checkpoint arrays use gzip/JSON with explicit dtype and shape. Undefined ancillary diagnostics, such as the boundary residual of a constant field, use typed nonfinite-float markers in otherwise strict JSON and recover their numerical NaN/Infinity values on reading. No executable pickle serialization is used.

The bundled historical checkpoints retain the SHA-256 of an acquired source copy with one extra final newline. The repository's core scripts match their authoritative GitHub blobs. Numerical verification checks source identity separately and excludes historical provenance hashes and elapsed times from scientific equality comparisons. New checkpoints record the current source hash.
