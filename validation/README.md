# Added uncertainty and sensitivity analyses

These wrappers import the original root-level detector and REA functions. They stream exact integer-level image means and evaluate spectra from one image at a time; the scientific algorithms are unchanged. See the repository README for actual parameters, results and limitations.

Run from the repository root:

```bash
python validation/check_reference_results.py
python validation/rea_validation_wrapper.py --aggregate --bootstrap 10000
python validation/rstat_sensitivity_wrapper.py --aggregate
python validation/rebuild_derived.py
```

All aggregation commands default to `validation/checkpoints/` and write to `outputs/validation/`. They work without raw images. The first command runs the three others into a temporary directory and verifies every scientific value against the bundled reference results, including all bootstrap draw indices and selections. Reference and sensitivity JSON files use strict JSON null for undefined/non-selected results. In the bootstrap NPZ only, `selected_pixels=inf` marks non-selection.

To recompute from original data, use `--image 1` through `--image 7` for each wrapper, with `REA_DATA_ROOT` or `--data-root` pointing to a folder containing `BSE/` and `reference_outputs/`. Always pass `--checkpoint-dir outputs/validation/checkpoints` to the image calls **and** each subsequent aggregation call. Image calls otherwise write output-dir/checkpoints while aggregation otherwise reads the bundled reference checkpoints. `--output-dir` is independent of the checkpoint location. `rebuild_derived.py` also accepts `--checkpoint-dir`.

`checkpoints/image_*.json.gz` contain baseline, masked49 and unmasked121 spectra, sample means, candidate positions, masks' coverage and availability. `checkpoints/rstat_image_*.json.gz` contain the 1536/3072 screening-window experiments. They were losslessly converted from trusted original local checkpoints; `checkpoint_io.py` uses gzip/JSON with explicit array dtype/shape and never executes pickle objects. Newly generated checkpoints use the same format.

The source SHA-256 in historical checkpoints/results belongs to the acquired source copy with one extra final newline. The released root scripts match their original GitHub blob hashes exactly. The verification excludes provenance hashes and elapsed times from numerical equality checks and independently checks the authoritative Git blob IDs. Both versions execute the same code.

`results/reference_sensitivity_extra.json` includes the shorter reference checks and a separately flagged exploratory 2560-pixel result that bypasses the original 49-position availability gate. That exploratory branch is not a production-pipeline result and is not a manuscript conclusion. `rebuild_derived.py` recreates that file and the Figure 2 input NPZ using the original ensemble-spectrum routine.

`inputs/original_combined_window_statistics.csv` is an unchanged copy of the structural CSV in the original Zenodo reference archive. It permits the aggregate baseline comparison without downloading the large images. Original images are required for recomputation of checkpoints; original numerical QEMSCAN maps are required for rerunning FEM. The added package does not claim a new FEM run.
