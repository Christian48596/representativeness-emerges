# Representativeness emerges from statistical admissibility and long-wavelength structural convergence

[![Dataset DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21717612.svg)](https://doi.org/10.5281/zenodo.21717612)

This repository contains the Python workflow used to:

- detect statistically admissible domains in two-dimensional material images;
- determine a structural representative elementary area (REA) from full-resolution BSE images;
- transfer the structurally predicted support to QEMSCAN-derived property maps;
- validate that support through periodic thermal and elastic homogenization.

The numerical validation reported in the associated manuscript is two-dimensional, so **REA** is used throughout the software and documentation. The corresponding three-dimensional concept is the representative elementary volume (REV).

## Data availability

The complete BSE and QEMSCAN datasets are archived separately on Zenodo:

**Dataset DOI:** [10.5281/zenodo.21717612](https://doi.org/10.5281/zenodo.21717612)

The Zenodo record contains the full-resolution BSE images, numerical QEMSCAN-derived property maps, stationary-domain masks, reference outputs, run logs, checksums, and detailed data documentation.

Large scientific data files are intentionally not stored in this GitHub repository.

## Repository contents

```text
representativeness-emerges/
├── README.md
├── README_DATA.md
├── requirements.txt
├── .gitignore
├── detect_stationary_domains_general.py
├── rea_stationary_mask_cli_v24_general.py
└── property_full_homogenized_inside_property_masks.py
```

The three executable scripts are:

| Script | Purpose |
|---|---|
| `detect_stationary_domains_general.py` | Detects statistically admissible domains in BSE images or numerical property maps. |
| `rea_stationary_mask_cli_v24_general.py` | Determines the structural REA from BSE images restricted to stationary-domain masks. |
| `property_full_homogenized_inside_property_masks.py` | Computes apparent thermal and elastic properties inside stationary QEMSCAN property masks. |

`README_DATA.md` describes the Zenodo archives, input units, expected directory structure, and data-specific limitations.

## Requirements

Python 3.10 or newer is recommended.

Install the required packages with:

```bash
python3 -m pip install -r requirements.txt
```

A clean virtual environment can be created with:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```

Check the command-line interfaces with:

```bash
python3 detect_stationary_domains_general.py --help
python3 rea_stationary_mask_cli_v24_general.py --help
python3 property_full_homogenized_inside_property_masks.py --help
```

## Download and arrange the data

Download the dataset from the Zenodo DOI above.

After extracting the Zenodo archives, use the following working structure:

```text
representativeness-emerges/
├── detect_stationary_domains_general.py
├── rea_stationary_mask_cli_v24_general.py
├── property_full_homogenized_inside_property_masks.py
├── data/
│   ├── BSE/
│   │   ├── 1_cut.tif
│   │   ├── 2_cut.tif
│   │   ├── 3_cut.tif
│   │   ├── 4_cut.tif
│   │   ├── 5_cut.tif
│   │   ├── 6_cut.tif
│   │   └── 7_cut.tif
│   └── QEMSCAN/
│       ├── kappaMap_1.csv
│       ├── ...
│       ├── kappaMap_7.csv
│       ├── EMap_1.csv
│       ├── ...
│       ├── EMap_7.csv
│       ├── nuMap_1.csv
│       ├── ...
│       └── nuMap_7.csv
└── outputs/
```

The QEMSCAN CSV files must contain the numerical physical values:

```text
kappaMap_i.csv   thermal conductivity, κ(x,y), in W m^-1 K^-1
EMap_i.csv       Young's modulus, E(x,y), in GPa
nuMap_i.csv      Poisson ratio, ν(x,y), dimensionless
```

Do not replace the numerical CSV maps with TIFF or PNG visualization images produced through intensity rescaling. Such images do not preserve the physical values required by the homogenization solvers.

## Reproducing the analysis

Run the four steps in the order given below.

---

## 1. Detect statistically admissible BSE domains

This step scans each full-resolution BSE image using local square windows. A window is accepted when its local mean and standard deviation are compatible with robust image-level reference values.

```bash
nohup python3 detect_stationary_domains_general.py \
  data/BSE/1_cut.tif \
  data/BSE/2_cut.tif \
  data/BSE/3_cut.tif \
  data/BSE/4_cut.tif \
  data/BSE/5_cut.tif \
  data/BSE/6_cut.tif \
  data/BSE/7_cut.tif \
  --input-type image \
  --field-transform raw \
  --output-dir outputs/stationary_domains_FULLRES_R2048_STRICT \
  --max-side 0 \
  --stationarity-window 2048 \
  --stride 512 \
  --tau-mu 0.02 \
  --tau-sigma 0.06 \
  --min-component-windows 4 \
  --dpi 300 \
  > outputs/run_detect_stationary_domains_FULLRES_R2048_STRICT.log 2>&1 &
```

Important settings:

- `--max-side 0` preserves the original BSE resolution.
- `--stationarity-window 2048` uses local windows of \(2048 \times 2048\) pixels.
- `--stride 512` moves the local window by 512 pixels.
- `--tau-mu 0.02` sets the local-mean compatibility tolerance.
- `--tau-sigma 0.06` sets the local-standard-deviation compatibility tolerance.
- `--min-component-windows 4` removes small isolated accepted components.

The masks used by the next step follow the naming convention:

```text
1_cut_stationary_mask.tif
2_cut_stationary_mask.tif
...
7_cut_stationary_mask.tif
```

Monitor the run with:

```bash
tail -f outputs/run_detect_stationary_domains_FULLRES_R2048_STRICT.log
```

---

## 2. Determine the structural REA from the BSE images

This step treats the seven BSE images as independent fields of view. Candidate windows are retained only when at least 98% of their pixels lie inside the stationary-domain mask produced in Step 1.

```bash
nohup python3 rea_stationary_mask_cli_v24_general.py \
  data/BSE/1_cut.tif \
  data/BSE/2_cut.tif \
  data/BSE/3_cut.tif \
  data/BSE/4_cut.tif \
  data/BSE/5_cut.tif \
  data/BSE/6_cut.tif \
  data/BSE/7_cut.tif \
  --output-dir outputs/combined_REA_final_STRICT_v24 \
  --output-prefix figure4_REA_stationary_domains_STRICT_v24 \
  --representative-image 1 \
  --max-side 0 \
  --stationary-mask-dir outputs/stationary_domains_FULLRES_R2048_STRICT \
  --stationary-mask-suffix _stationary_mask.tif \
  --min-mask-fraction 0.98 \
  --window-sizes 256,384,512,640,768,896,1024,1280,1536,1792,2048 \
  --reference-L 2048 \
  --selected-spectra 256,512,1024,2048 \
  --n-common-k-grid 600 \
  --field-mode raw \
  --n-window-samples 121 \
  --errorbar-mode sem \
  --spectral-band low \
  --low-k-fraction 0.25 \
  --low-k-weight inverse \
  --tau-Z 0.025 \
  --tau-C 0.30 \
  --tau-ens 0.10 \
  --panel-b-mode raw+smooth \
  --panel-b-smooth-window 3 \
  --panel-c-mode raw+smooth \
  --panel-c-smooth-window 3 \
  --dpi 600 \
  --no-pdf \
  > outputs/run_REA_stationary_STRICT_v24.log 2>&1 &
```

The structural REA is selected from persistent convergence of:

- the scale-dependent mean BSE contrast;
- the low-wavenumber spectral residual;
- reproducibility across the seven independent images.

The display smoothing options modify only the plotted curves. They do not change the numerical REA decision.

For the manuscript dataset:

```text
L_REA = 1536 BSE pixels ≈ 2.01 mm
```

Monitor the run with:

```bash
tail -f outputs/run_REA_stationary_STRICT_v24.log
```

---

## 3. Detect statistically admissible QEMSCAN property domains

This step constructs stationary masks directly on the QEMSCAN property-map grid. Thermal conductivity, Young's modulus, and Poisson ratio are combined into one standardized stationarity field.

```bash
nohup python3 detect_stationary_domains_general.py \
  --kappa-maps \
    data/QEMSCAN/kappaMap_1.csv \
    data/QEMSCAN/kappaMap_2.csv \
    data/QEMSCAN/kappaMap_3.csv \
    data/QEMSCAN/kappaMap_4.csv \
    data/QEMSCAN/kappaMap_5.csv \
    data/QEMSCAN/kappaMap_6.csv \
    data/QEMSCAN/kappaMap_7.csv \
  --E-maps \
    data/QEMSCAN/EMap_1.csv \
    data/QEMSCAN/EMap_2.csv \
    data/QEMSCAN/EMap_3.csv \
    data/QEMSCAN/EMap_4.csv \
    data/QEMSCAN/EMap_5.csv \
    data/QEMSCAN/EMap_6.csv \
    data/QEMSCAN/EMap_7.csv \
  --nu-maps \
    data/QEMSCAN/nuMap_1.csv \
    data/QEMSCAN/nuMap_2.csv \
    data/QEMSCAN/nuMap_3.csv \
    data/QEMSCAN/nuMap_4.csv \
    data/QEMSCAN/nuMap_5.csv \
    data/QEMSCAN/nuMap_6.csv \
    data/QEMSCAN/nuMap_7.csv \
  --stationarity-source combined \
  --combined-transform log-robust \
  --output-dir outputs/stationary_property_combined_STRICT \
  --stationarity-window 256 \
  --stride 64 \
  --tau-mu 0.02 \
  --tau-sigma 0.06 \
  --min-component-windows 4 \
  --dpi 300 \
  > outputs/run_detect_stationary_property_combined_STRICT.log 2>&1 &
```

The expected masks are:

```text
property_combined_1_stationary_mask.tif
property_combined_2_stationary_mask.tif
...
property_combined_7_stationary_mask.tif
```

Monitor the run with:

```bash
tail -f outputs/run_detect_stationary_property_combined_STRICT.log
```

---

## 4. Validate the transferred REA through property homogenization

This step computes apparent thermal conductivity, elastic stiffness, and directional Young's moduli inside windows accepted by the QEMSCAN stationary masks.

```bash
nohup python3 property_full_homogenized_inside_property_masks.py \
  --kappa-maps \
    data/QEMSCAN/kappaMap_1.csv \
    data/QEMSCAN/kappaMap_2.csv \
    data/QEMSCAN/kappaMap_3.csv \
    data/QEMSCAN/kappaMap_4.csv \
    data/QEMSCAN/kappaMap_5.csv \
    data/QEMSCAN/kappaMap_6.csv \
    data/QEMSCAN/kappaMap_7.csv \
  --E-maps \
    data/QEMSCAN/EMap_1.csv \
    data/QEMSCAN/EMap_2.csv \
    data/QEMSCAN/EMap_3.csv \
    data/QEMSCAN/EMap_4.csv \
    data/QEMSCAN/EMap_5.csv \
    data/QEMSCAN/EMap_6.csv \
    data/QEMSCAN/EMap_7.csv \
  --nu-maps \
    data/QEMSCAN/nuMap_1.csv \
    data/QEMSCAN/nuMap_2.csv \
    data/QEMSCAN/nuMap_3.csv \
    data/QEMSCAN/nuMap_4.csv \
    data/QEMSCAN/nuMap_5.csv \
    data/QEMSCAN/nuMap_6.csv \
    data/QEMSCAN/nuMap_7.csv \
  --mask-dir outputs/stationary_property_combined_STRICT \
  --mask-template 'property_combined_{i}_stationary_mask.tif' \
  --output-dir outputs/property_full_homogenized_property_masks_STRICT \
  --output-prefix figure5_property_full_homogenized_property_masks_STRICT \
  --window-sizes 32,48,64,80,96,112,128,160,192,204,224,256,320,384,448,512 \
  --L-rea 204 \
  --property-pixel-size-mm 0.01 \
  --min-mask-fraction 0.98 \
  --n-window-samples 4 \
  --elastic-mode plane_strain \
  --errorbar-mode sem \
  --dpi 600 \
  --no-pdf \
  > outputs/run_property_full_homogenized_property_masks_STRICT.log 2>&1 &
```

The option:

```text
--L-rea 204
```

marks the structural REA transferred from the BSE grid to the QEMSCAN property grid. It is an independently predicted support and is not fitted to the apparent thermal or elastic property curves.

With a QEMSCAN property-map pixel size of \(0.01~\mathrm{mm/pixel}\):

```text
204 pixels ≈ 2.04 mm
```

The property script writes:

```text
*_kappa_app_curves.png
*_C_app_curves.png
*_E_app_curves.png
*_summary.csv
*_window_values.csv
```

Monitor the run with:

```bash
tail -f outputs/run_property_full_homogenized_property_masks_STRICT.log
```

## Workflow summary

```text
BSE images
   |
   v
stationary-domain detection
   |
   v
BSE stationary masks
   |
   v
structural REA analysis
   |
   v
L_REA = 1536 BSE pixels ≈ 2.01 mm
   |
   | transfer by physical length
   v
L_REA = 204 QEMSCAN pixels ≈ 2.04 mm
   |
   v
QEMSCAN stationary masks
   |
   v
periodic thermal and elastic homogenization
   |
   v
independent property-level validation
```

## Reproducibility notes

- Keep the ordering of the seven BSE images consistent throughout the workflow.
- Keep the ordering of `kappaMap_i.csv`, `EMap_i.csv`, and `nuMap_i.csv` identical.
- Property maps and their corresponding masks must have exactly the same dimensions.
- Use the numerical QEMSCAN CSV maps in physical units.
- `--max-side 0` is required to reproduce the full-resolution BSE analysis.
- Candidate windows are accepted only when their stationary-mask coverage satisfies `--min-mask-fraction`.
- The smoothed curves are display aids and do not enter the structural REA criterion.
- The property-level elastic calculation is the most computationally expensive stage.
- The property script uses a deterministic random seed by default unless `--rng-seed` is changed.
- Small floating-point differences may occur across operating systems, BLAS libraries, and sparse linear solvers.

## Reference outputs

Reference masks, figures, numerical summaries, and run logs are provided in the Zenodo dataset. These archived outputs are intended for comparison with newly generated results.

## Citation

When using this repository, cite:

1. the associated manuscript;
2. the Zenodo dataset using the citation displayed on the Zenodo record;
3. the software release DOI, once a GitHub release has been archived separately on Zenodo.

Dataset:

```text
https://doi.org/10.5281/zenodo.21717612
```

## License

The dataset license is specified in the Zenodo record.

The software license should be stated in a separate `LICENSE` file in this repository. Confirm the selected software license with all code authors before publishing the repository.

## Contact

```text
Dr. Christian Tantardini
christiantantardini@ymail.com

Prof. Eduardo Garzanti
eduardo.garzanti@unimib.it
```
