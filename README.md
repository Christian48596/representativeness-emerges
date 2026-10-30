# Statistical admissibility and long-wavelength structural convergence determine representative support in granular materials

This repository contains the computational workflow and reproducibility files for the manuscript by Christian Tantardini, Fernando Alonso-Marroquín, Abdullah Alqubalee and Eduardo Garzanti. The demonstrated system is dune sand imaged in two dimensions. The reported representative elementary area (REA) is a **conditional, finite-reference structural support estimate**, not a global representative volume or a universal material constant.

The three original computational scripts are unchanged from commit [`8bd49ef71bbd8e320bdb00b3633f8d137cb72bb0`](https://github.com/Christian48596/representativeness-emerges/tree/8bd49ef71bbd8e320bdb00b3633f8d137cb72bb0). The additions provide uncertainty and sensitivity analyses, corrected documentation, compact numerical inputs, and publication figure layouts. No original selection or homogenization algorithm has been changed.

## Data and contents

Full-resolution BSE images, numerical QEMSCAN property maps, masks and reference outputs are available from [Zenodo concept DOI 10.5281/zenodo.21717612](https://doi.org/10.5281/zenodo.21717612). The original inputs used here are version [10.5281/zenodo.21717613](https://doi.org/10.5281/zenodo.21717613). The updated `reference_outputs_and_masks.tar.gz` contains the new derived results and publication figures alongside the original outputs. Cite the specific Zenodo version actually downloaded and the Git commit/release actually used; a new version DOI is assigned by Zenodo when that version is published.

| Location | Contents |
|---|---|
| Three root-level Python scripts | Original detector, structural REA analysis and periodic property homogenization |
| `docs/ORIGINAL_RUN_COMMANDS.md` | Explicit original-run commands and implementation qualifications |
| `validation/` | External validation wrappers and compact verification command |
| `validation/checkpoints/` | Safe gzip/JSON per-image spectra and statistics; no executable pickle files |
| `validation/results/` | Baseline, bootstrap, sensitivity results and seven per-image CSV tables |
| `validation/inputs/` | Original archived structural CSV used to check the baseline |
| `figures/` | Scripts for revised Figures 1–3 and reconstructed QEMSCAN maps |
| `derived_data/` | Compact verified CSV/NPZ inputs for figure regeneration |

Raw images and full-resolution masks belong in the Zenodo data directory; they are not duplicated in GitHub. Extract the three original data archives into one directory containing `BSE/`, `QEMSCAN/` and `reference_outputs/`.

## Install and verify compact results

Use Python 3.10 or newer and the four dependencies in the original `requirements.txt`:

```bash
python -m pip install -r requirements.txt
python validation/check_reference_results.py
```

The verification command requires no raw images. It checks the Git blob hashes of all three original scripts, reruns all 10,000 bootstrap resamples from the bundled checkpoints, and compares the numerical baseline, leave-one-image-out results, sensitivity cases, seven per-image CSVs and Figure 2 inputs. Recorded wall times vary between machines and are excluded from equality checks. `requirements-revision.txt` provides the exact package versions used for release verification; these are not claimed to be the unknown environment of the original FEM run.

## Reproduce the added analyses

Run from the repository root. These commands use the bundled safe checkpoints and write new files under `outputs/validation/`:

```bash
python validation/rea_validation_wrapper.py --aggregate --bootstrap 10000
python validation/rstat_sensitivity_wrapper.py --aggregate
python validation/rebuild_derived.py
```

The first command calls the original aggregation and selection functions for the baseline, masked 49-position and unmasked 121-position comparisons, larger-reference availability checks, whole-image bootstrap and leave-one-image-out analysis. The second aggregates the alternative screening-window results. The third regenerates shorter-reference cases, explicitly identified exploratory cases and the Figure 2 numeric input archive.

To recompute the underlying statistics and alternative masks from the raw data, set the data root and process the seven images sequentially. Each image uses several gigabytes of memory because the raw field contains approximately 205 million pixels; avoid running seven image jobs concurrently on a small machine.

```bash
export REA_DATA_ROOT="/absolute/path/to/extracted/data"
for i in 1 2 3 4 5 6 7; do
  python validation/rea_validation_wrapper.py --image "$i" --checkpoint-dir outputs/validation/checkpoints
  python validation/rstat_sensitivity_wrapper.py --image "$i" --checkpoint-dir outputs/validation/checkpoints
done
python validation/rea_validation_wrapper.py --aggregate --bootstrap 10000 --checkpoint-dir outputs/validation/checkpoints
python validation/rstat_sensitivity_wrapper.py --aggregate --checkpoint-dir outputs/validation/checkpoints
python validation/rebuild_derived.py --checkpoint-dir outputs/validation/checkpoints
```

The explicit `--checkpoint-dir` in all three aggregation commands selects the freshly recomputed values. Without it, aggregation reads the bundled checkpoints. `--data-root` overrides the environment variable, and `--output-dir` changes the output destination. Original computational scripts do not read `REA_DATA_ROOT`; their explicit commands are in `docs/ORIGINAL_RUN_COMMANDS.md`.

## Regenerate figures

```bash
python figures/regenerate_figure1_layout.py --data-root "/absolute/path/to/extracted/data"
python figures/regenerate_figure2_layout.py --data-root "/absolute/path/to/extracted/data"
python figures/regenerate_revised_figures.py
```

Figures 1 and 2 require the original BSE fields (Figure 1 also reads its original mask and diagnostic CSV). Figure 3 and reconstructed QEMSCAN maps can be regenerated entirely from the compact inputs. See `figures/README.md` for output paths and the optional full-map route. The Figure 1/2 changes remove titles, use bold external panel letters and equal panel dimensions. Numerical spectra and convergence curves are preserved. Figure 3 reports linear-scale mineral area fractions and adds shear compliance information from inversion of each full archived elastic stiffness matrix.

## Actual settings and interpretation

The full-resolution BSE run uses support sides 256, 384, 512, 640, 768, 896, 1024, 1280, 1536, 1792 and 2048 pixels; reference side 2048; 121 nominal sampling positions; and 600 common radial wavenumbers. BSE admissibility uses a 2048-pixel screening window, stride 512, mean tolerance 0.02, standard-deviation tolerance 0.06 and minimum four coarse connected cells. These are explicit run parameters. CLI defaults (including downsampling, reference 896, 49 samples and 300 wavenumbers) do not reproduce the reported full-resolution run.

A separate preliminary availability gate in the original code uses 49 nominal positions even when actual sampling requests 121. Accepted analysis windows have at least 98% mask coverage. Complete rectangular windows are analyzed, including their permitted outside-mask pixels; the mask is not multiplied into the image. Local marginal mean/variance admissibility does not establish full spatial stationarity.

The executable selection uses the regularized mean residual and the unsmoothed spectral residual. Five-point median regularization and the suffix-maximum mean envelope enter selection. Optional three-point plotted trends are display aids. The auxiliary ensemble variability diagnostic does **not** enter acceptance. Persistent acceptance requires at least two nonreference supports and the implemented suffix-dependent tail count. Consequently, 1536 pixels is the largest selectable side in the reported sequence.

| Analysis | Result |
|---|---|
| Baseline | 1536 BSE pixels = 2.01216 mm, relative to reference 2048 pixels |
| 10,000 whole-image bootstrap resamples, seed 20261029 | 3310 non-selections (33.10%) |
| Successful bootstrap selections only | 2.5–97.5 percentile range 896–1536 pixels; **not an unconditional 95% confidence interval** |
| Masked 49-position sampling | No selected support |
| Unmasked 121-position sampling | 1280 pixels; this is a different sampling population |
| Screening window 1536, stride 512 | Only two common valid supports; no estimate |
| Screening window 3072, stride 512 | 1536-pixel estimate |
| Truncated reference 1536 / 1792 | No estimate / 1280 pixels |
| Reference 2560–4096 with production availability gate | At least one image fails eligibility |

Bootstrap resampling treats the seven images as independent clusters and holds the masks, support grid and reference fixed. It does not include uncertainty from threshold choice, mask estimation, resolution or segmentation. The finite reference is only 1.33 times the selected side and is not established as an asymptotic ground truth. Pixel size and microstructural resolution may change the result. The demonstrated domain and material do not establish general validity in three dimensions or other materials.

The file `reference_sensitivity_extra.json` separately retains an exploratory 2560-pixel-reference calculation that bypasses the original 49-position gate. Its `enforces_original49_common_size_gate` flag is false; it is **not** a result of the unmodified production eligibility pipeline and was not used as a manuscript conclusion.

Property calculations use 2D plane strain and in-plane conduction, periodic Q4 elements and direct sparse solves. Property masks retain images 1, 5 and 6; supports 32–320 contain four windows from each image, while 384 contains four windows from image 6 only. The marked 204-pixel property support is 2.04 mm. Pooled-window SEM is descriptive and is not an independent-image confidence interval. Full stiffness inversion retains normal/shear coupling; `Gxy=1/S33` differs from `C33`.

## Formats and provenance

JSON `null` denotes an undefined statistic, unavailable result or non-selection, as described by its field. In `bootstrap_draws.npz`, infinity in `selected_pixels` is an internal non-selection marker, not an infinite physical support. Load NPZ files with `allow_pickle=False`. Gzip/JSON checkpoint arrays carry dtype and shape and are loaded without executable serialization.

Historical checkpoints and result files retain the SHA-256 of the acquisition copy of the structural script. That copy differed from the GitHub file only by one extra final newline. This release contains the exact GitHub blobs; new aggregations record their current SHA-256. Numerical comparisons exclude only provenance hashes and elapsed times. The source verification command checks the three authoritative Git blob IDs.

This package is prepared for the manuscript revision; a suggested release label is `v1.1.0`. Preparation does not imply that a GitHub release or new Zenodo DOI already exists.

## License and contact

The repository contains the original GNU General Public License v3 text in `LICENSE`. Dataset licensing is specified by the Zenodo record. Contact: Dr. Christian Tantardini, christiantantardini@ymail.com; Prof. Eduardo Garzanti, eduardo.garzanti@unimib.it.
