# Original calculation commands

These commands document the archived calculation using the three unchanged computational scripts. They make the run parameters explicit and write fresh outputs under `outputs/original_run/`; the archived data and outputs remain intact. They are reproduction instructions, not a claim that the entire original workflow was rerun when this upload package was prepared.

Run the blocks in order in **Bash**, from the repository root, after installing `requirements.txt`. Extract the three Zenodo archives into one directory containing `BSE/`, `QEMSCAN/` and `reference_outputs/`. Set its absolute path below. Commands run in the foreground. Check that each command succeeds before starting the next stage.

```bash
export REA_DATA_ROOT="/absolute/path/to/extracted/data"
REA_RUN_DIR="$PWD/outputs/original_run"
mkdir -p "$REA_RUN_DIR"

bse_files=()
kappa_files=()
elastic_files=()
poisson_files=()
for i in {1..7}; do
  bse_files+=("$REA_DATA_ROOT/BSE/${i}_cut.tif")
  kappa_files+=("$REA_DATA_ROOT/QEMSCAN/kappaMap_${i}.csv")
  elastic_files+=("$REA_DATA_ROOT/QEMSCAN/EMap_${i}.csv")
  poisson_files+=("$REA_DATA_ROOT/QEMSCAN/nuMap_${i}.csv")
done
```

## 1. BSE admissibility masks

```bash
python3 detect_stationary_domains_general.py "${bse_files[@]}" \
  --input-type image --field-transform raw --max-side 0 \
  --output-dir "$REA_RUN_DIR/stationary_domains_FULLRES_R2048_STRICT" \
  --stationarity-window 2048 --stride 512 \
  --tau-mu 0.02 --tau-sigma 0.06 --min-component-windows 4 --dpi 300 \
  > "$REA_RUN_DIR/run_detect_stationary_domains_FULLRES_R2048_STRICT.log" 2>&1
```

The source's morphology operations remain enabled. The output masks are named `1_cut_stationary_mask.tif` through `7_cut_stationary_mask.tif`. This is a local marginal-compatibility screen, not a proof of full spatial stationarity.

## 2. Structural support selection

```bash
python3 rea_stationary_mask_cli_v24_general.py "${bse_files[@]}" \
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

The archived result is 1536 pixels (2.01216 mm using 1.31 micrometres per BSE pixel). The reference is **2048 pixels**, not 4096. `--pixel-size 1` preserves the actual pixel-coordinate spectral convention; physical reporting of lengths is a separate conversion. The requested 121 sampling positions form an 11-by-11 grid, plus the centre if distinct. A separate source check of common-support availability uses 49 positions. Ensemble reproducibility (`--tau-ens`, also accepted under the legacy name `--tau-boundary`) is an auxiliary diagnostic, not another selection gate. Plot smoothing differs from the source's five-point median and suffix-maximum regularization used in selection.

Use a basename for `--output-prefix`, as above: the script joins it to `--output-dir`.

## 3. QEMSCAN admissibility masks

```bash
python3 detect_stationary_domains_general.py \
  --kappa-maps "${kappa_files[@]}" \
  --E-maps "${elastic_files[@]}" \
  --nu-maps "${poisson_files[@]}" \
  --stationarity-source combined --combined-transform log-robust --max-side 0 \
  --output-dir "$REA_RUN_DIR/stationary_property_combined_STRICT" \
  --stationarity-window 256 --stride 64 \
  --tau-mu 0.02 --tau-sigma 0.06 --min-component-windows 4 --dpi 300 \
  > "$REA_RUN_DIR/run_detect_stationary_property_combined_STRICT.log" 2>&1
```

Use the numerical CSV property maps. Intensity-rescaled image exports do not preserve physical property values. The combined-field screen produces `property_combined_1_stationary_mask.tif` through `property_combined_7_stationary_mask.tif`.

## 4. Periodic property homogenization

```bash
python3 property_full_homogenized_inside_property_masks.py \
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

The seed is the unchanged source default, made explicit. The nominal transferred support is 204 pixels (2.04 mm). The archived calculation contains 160 solved windows: 12 per support from 32 to 320 pixels, drawn from images 1, 5 and 6, and four at 384 pixels from image 6 only. No accepted windows are available at 448 or 512 pixels. Elasticity is two-dimensional plane strain; conduction is in-plane. Pooled-window SEM is descriptive, not an independent-image confidence interval.

The original programs generate their historical plot layouts. For current publication figures, full-compliance shear postprocessing, bootstrap uncertainty and sensitivity analyses, use the added commands in the [repository README](../README.md), [validation documentation](../validation/README.md) and [figure documentation](../figures/README.md). These added analyses do not change the three original computational scripts.
