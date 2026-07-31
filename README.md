# Representativeness emerges from statistical admissibility and long-wavelength structural convergence

This repository contains the Python workflow used to identify statistically admissible domains, determine a structural representative elementary area (REA) from BSE images, and validate the predicted support using thermal and elastic property fields derived from QEMSCAN mineral maps.

The workflow has two connected parts:

1. **BSE structural analysis**
   - detect statistically compatible regions in each BSE image;
   - restrict finite-window sampling to those regions;
   - determine the structural REA from persistent convergence of the image mean and the low-wavenumber covariance spectrum.

2. **QEMSCAN property validation**
   - construct stationary masks on the property-map grid from thermal and elastic fields;
   - sample only windows that lie almost entirely inside those masks;
   - compute apparent conductivity, stiffness, and directional Young's moduli using periodic finite-element homogenization;
   - compare their size dependence with the structural support predicted from BSE images.

The numerical validation in the associated manuscript is two-dimensional. The general theoretical framework is formulated for both two-dimensional images and three-dimensional voxelized volumes.

---

### Repository contents

```text
.
├── detect_stationary_domains_general.py
├── rev_stationary_mask_cli_v24_general.py
├── property_full_homogenized_inside_property_masks.py
├── README.md
├── BSE/
│   ├── 1_cut.tif
│   ├── 2_cut.tif
│   ├── ...
│   └── 7_cut.tif
└── QEMSCAN/
    ├── kappaMap_1.csv
    ├── ...
    ├── kappaMap_7.csv
    ├── EMap_1.csv
    ├── ...
    ├── EMap_7.csv
    ├── nuMap_1.csv
    ├── ...
    └── nuMap_7.csv
```

The commands below assume that the scripts and input files are in the current working directory. When the data are stored in subdirectories, replace each filename with its relative or absolute path.

The Python files should use these exact names:

```text
detect_stationary_domains_general.py
rev_stationary_mask_cli_v24_general.py
property_full_homogenized_inside_property_masks.py
```

---

### Requirements

A recent Python 3 installation is required. Python 3.10 or newer is recommended.

Install the required packages with:

```bash
python -m pip install numpy scipy matplotlib pillow
```

A clean virtual environment can be created with:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy scipy matplotlib pillow
```

On Windows PowerShell, activate the environment with:

```powershell
.venv\Scripts\Activate.ps1
```

Check that each program is available:

```bash
python detect_stationary_domains_general.py --help
python rev_stationary_mask_cli_v24_general.py --help
python property_full_homogenized_inside_property_masks.py --help
```

---

### Input data

#### BSE images

The BSE analysis uses seven grayscale TIFF images:

```text
1_cut.tif
2_cut.tif
3_cut.tif
4_cut.tif
5_cut.tif
6_cut.tif
7_cut.tif
```

Gray intensity is treated as a local compositional or average-atomic-number proxy. It is not interpreted as thermal conductivity or as a direct quantitative chemical concentration.

The manuscript calculations use the full image resolution. This is selected with:

```text
--max-side 0
```

#### QEMSCAN-derived property maps

For every mapped section, three numerical arrays are required:

```text
kappaMap_i.csv   thermal conductivity, W m^-1 K^-1
EMap_i.csv       Young's modulus, GPa
nuMap_i.csv      Poisson ratio, dimensionless
```

The files with the same index must describe the same section and must have identical array dimensions. For example:

```text
kappaMap_1.csv
EMap_1.csv
nuMap_1.csv
```

must form one aligned property set.

Use the **numerical CSV maps containing physical values**. Do not use TIFF or PNG visualization maps created through intensity rescaling, such as `mat2gray` or integer image conversion, because those files do not preserve the physical units required for homogenization.

---

### Complete workflow

Run the four steps in the order shown below.

---

### Step 1: Detect stationary domains in the BSE images

**Script:** `detect_stationary_domains_general.py`

This step scans each full-resolution BSE image using local square windows. For every window, it calculates the local mean and standard deviation and compares them with robust image-level reference values. Windows satisfying both tolerances are converted into a binary stationary-domain mask.

**Run:**

```bash
nohup python detect_stationary_domains_general.py \
  1_cut.tif 2_cut.tif 3_cut.tif 4_cut.tif 5_cut.tif 6_cut.tif 7_cut.tif \
  --input-type image \
  --field-transform raw \
  --output-dir stationary_domains_FULLRES_R2048_STRICT \
  --max-side 0 \
  --stationarity-window 2048 \
  --stride 512 \
  --tau-mu 0.02 \
  --tau-sigma 0.06 \
  --min-component-windows 4 \
  --dpi 300 > run_detect_stationary_domains_FULLRES_R2048_STRICT.log 2>&1 &
```

**Main settings:**

- `--max-side 0` keeps the original image resolution.
- `--stationarity-window 2048` uses local windows of \(2048\times2048\) pixels.
- `--stride 512` moves the stationarity window by 512 pixels.
- `--tau-mu 0.02` allows a relative local-mean deviation of 2%.
- `--tau-sigma 0.06` allows a relative local-standard-deviation deviation of 6%.
- `--min-component-windows 4` removes isolated accepted regions smaller than four coarse stationarity windows.
- `--field-transform raw` uses the original BSE gray levels without logarithmic or standardized transformation.

**Output directory:**

```text
stationary_domains_FULLRES_R2048_STRICT/
```

For each BSE image, the script writes files such as:

```text
1_cut_stationary_mask.tif
1_cut_stationary_mask.png
1_cut_stationary_overlay.png
1_cut_stationary_summary.png
1_cut_stationary_windows.csv
```

It also writes:

```text
stationary_domain_global_summary.csv
```

The TIFF masks are used directly in Step 2.

**Monitor the run:**

```bash
tail -f run_detect_stationary_domains_FULLRES_R2048_STRICT.log
```

---

### Step 2: Determine the structural REA from the BSE images

**Script:** `rev_stationary_mask_cli_v24_general.py`

This step treats the seven BSE images as independent realizations of the same regional material population. Candidate square windows are accepted only when at least 98% of their pixels lie inside the stationary-domain mask produced in Step 1.

The structural support is selected from persistent convergence of:

- the apparent gray-level mean;
- the low-wavenumber covariance-spectrum residual.

The low-wavenumber region is emphasized because it contains the longest spatial fluctuations and normally converges more slowly than local image features.

**Run:**

```bash
nohup python rev_stationary_mask_cli_v24_general.py \
  1_cut.tif 2_cut.tif 3_cut.tif 4_cut.tif 5_cut.tif 6_cut.tif 7_cut.tif \
  --output-dir combined_REV_final_STRICT_v24 \
  --output-prefix figure4_REV_stationary_domains_STRICT_v24 \
  --representative-image 1 \
  --max-side 0 \
  --stationary-mask-dir stationary_domains_FULLRES_R2048_STRICT \
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
  --tau-boundary 0.10 \
  --panel-b-mode raw+smooth \
  --panel-b-smooth-window 3 \
  --panel-c-mode raw+smooth \
  --panel-c-smooth-window 3 \
  --dpi 600 \
  --no-pdf > run_REV_stationary_STRICT_v24.log 2>&1 &
```

**Main settings:**

- `--stationary-mask-dir` points to the masks generated in Step 1.
- `--stationary-mask-suffix _stationary_mask.tif` matches masks to images by filename stem.
- `--min-mask-fraction 0.98` requires at least 98% stationary-domain coverage.
- `--window-sizes` gives all tested BSE support sizes in pixels.
- `--reference-L 2048` defines the largest finite reference support.
- `--selected-spectra 256,512,1024,2048` selects the spectra displayed in the summary figure.
- `--n-window-samples 121` requests broad sampling of admissible windows.
- `--errorbar-mode sem` reports standard errors of the mean.
- `--spectral-band low` restricts the structural comparison to low wavenumbers.
- `--low-k-fraction 0.25` uses the lowest 25% of the common radial wavenumber grid.
- `--low-k-weight inverse` gives greater weight to the longest resolved spatial fluctuations.
- `--tau-Z 0.025` sets the apparent-mean tolerance.
- `--tau-C 0.30` sets the covariance-spectrum tolerance.
- `--tau-boundary 0.10` controls the auxiliary boundary-consistency diagnostic.
- `--panel-b-mode raw+smooth` and `--panel-c-mode raw+smooth` add running-median trends for visual clarity only. The numerical REA decision remains based on the unsmoothed values.
- `--no-pdf` saves PNG output only.

**Output directory:**

```text
combined_REV_final_STRICT_v24/
```

The output includes the manuscript-style four-panel structural summary, a separate consistency-check figure, numerical convergence tables, and run diagnostics. The exact filenames begin with:

```text
figure4_REV_stationary_domains_STRICT_v24
```

For the manuscript dataset, the selected structural support is:

```text
L_REA = 1536 BSE pixels
```

Using the BSE pixel size of approximately \(1.31~\mu\mathrm{m}\), this corresponds to a physical side length of approximately:

```text
2.01 mm
```

**Monitor the run:**

```bash
tail -f run_REV_stationary_STRICT_v24.log
```

---

### Step 3: Detect stationary domains on the QEMSCAN property grid

**Script:** `detect_stationary_domains_general.py`

This step builds stationary masks directly on the QEMSCAN property-map grid. The conductivity, Young's-modulus, and Poisson-ratio channels are robustly standardized and combined into one scalar stationarity field.

Conductivity and Young's modulus are logarithmically transformed before robust standardization because their phase contrasts can span a wide range. Poisson's ratio is standardized without a logarithmic transform.

**Run:**

```bash
python detect_stationary_domains_general.py \
  --kappa-maps kappaMap_1.csv kappaMap_2.csv kappaMap_3.csv kappaMap_4.csv kappaMap_5.csv kappaMap_6.csv kappaMap_7.csv \
  --E-maps EMap_1.csv EMap_2.csv EMap_3.csv EMap_4.csv EMap_5.csv EMap_6.csv EMap_7.csv \
  --nu-maps nuMap_1.csv nuMap_2.csv nuMap_3.csv nuMap_4.csv nuMap_5.csv nuMap_6.csv nuMap_7.csv \
  --stationarity-source combined \
  --combined-transform log-robust \
  --output-dir stationary_property_combined_STRICT \
  --stationarity-window 256 \
  --stride 64 \
  --tau-mu 0.02 \
  --tau-sigma 0.06 \
  --min-component-windows 4 \
  --dpi 300
```

**Main settings:**

- `--stationarity-source combined` uses all three property channels.
- `--combined-transform log-robust` logarithmically transforms conductivity and Young's modulus, robustly standardizes all channels, and combines them into one stationarity field.
- `--stationarity-window 256` uses \(256\times256\)-pixel local windows on the property grid.
- `--stride 64` moves the stationarity window by 64 property-map pixels.
- `--tau-mu 0.02` and `--tau-sigma 0.06` apply the same strict local compatibility tolerances used for the BSE masks.
- `--min-component-windows 4` removes small isolated stationary components.

**Output directory:**

```text
stationary_property_combined_STRICT/
```

The masks required by Step 4 follow this naming pattern:

```text
property_combined_1_stationary_mask.tif
property_combined_2_stationary_mask.tif
...
property_combined_7_stationary_mask.tif
```

The masks remain on the same numerical grid as the QEMSCAN-derived property maps.

To run this step in the background, use:

```bash
nohup python detect_stationary_domains_general.py \
  --kappa-maps kappaMap_1.csv kappaMap_2.csv kappaMap_3.csv kappaMap_4.csv kappaMap_5.csv kappaMap_6.csv kappaMap_7.csv \
  --E-maps EMap_1.csv EMap_2.csv EMap_3.csv EMap_4.csv EMap_5.csv EMap_6.csv EMap_7.csv \
  --nu-maps nuMap_1.csv nuMap_2.csv nuMap_3.csv nuMap_4.csv nuMap_5.csv nuMap_6.csv nuMap_7.csv \
  --stationarity-source combined \
  --combined-transform log-robust \
  --output-dir stationary_property_combined_STRICT \
  --stationarity-window 256 \
  --stride 64 \
  --tau-mu 0.02 \
  --tau-sigma 0.06 \
  --min-component-windows 4 \
  --dpi 300 > run_detect_stationary_property_combined_STRICT.log 2>&1 &
```

Monitor it with:

```bash
tail -f run_detect_stationary_property_combined_STRICT.log
```

---

### Step 4: Validate the structural support by property homogenization

**Script:** `property_full_homogenized_inside_property_masks.py`

This step solves periodic finite-element cell problems inside windows accepted by the QEMSCAN property masks.

The script computes:

- apparent thermal-conductivity tensor components;
- apparent two-dimensional elastic stiffness components;
- directional apparent Young's moduli derived from the apparent compliance matrix.

The candidate windows are accepted only when at least 98% of their pixels lie inside the property-grid stationary mask.

**Run:**

```bash
nohup python property_full_homogenized_inside_property_masks.py \
  --kappa-maps kappaMap_1.csv kappaMap_2.csv kappaMap_3.csv kappaMap_4.csv kappaMap_5.csv kappaMap_6.csv kappaMap_7.csv \
  --E-maps EMap_1.csv EMap_2.csv EMap_3.csv EMap_4.csv EMap_5.csv EMap_6.csv EMap_7.csv \
  --nu-maps nuMap_1.csv nuMap_2.csv nuMap_3.csv nuMap_4.csv nuMap_5.csv nuMap_6.csv nuMap_7.csv \
  --mask-dir stationary_property_combined_STRICT \
  --mask-template 'property_combined_{i}_stationary_mask.tif' \
  --output-dir property_full_homogenized_property_masks_STRICT \
  --output-prefix figure5_property_full_homogenized_property_masks_STRICT \
  --window-sizes 32,48,64,80,96,112,128,160,192,204,224,256,320,384,448,512 \
  --L-rev 204 \
  --property-pixel-size-mm 0.01 \
  --min-mask-fraction 0.98 \
  --n-window-samples 4 \
  --elastic-mode plane_strain \
  --errorbar-mode sem \
  --dpi 600 \
  --no-pdf > run_property_full_homogenized_property_masks_STRICT.log 2>&1 &
```

**Main settings:**

- `--mask-dir` points to the property-grid masks generated in Step 3.
- `--mask-template 'property_combined_{i}_stationary_mask.tif'` maps property set \(i\) to its stationary mask.
- `--window-sizes` gives all tested support sizes in property-map pixels.
- `--L-rev 204` places the BSE-derived structural prediction on the property-grid convergence plots.
- `--property-pixel-size-mm 0.01` converts property-map pixels to millimetres.
- `--min-mask-fraction 0.98` requires at least 98% stationary-mask coverage.
- `--n-window-samples 4` uses four accepted windows per image and support size. The full elastic calculation is computationally expensive.
- `--elastic-mode plane_strain` uses the two-dimensional plane-strain constitutive assumption.
- `--errorbar-mode sem` reports standard errors across accepted windows and images.
- `--no-pdf` saves PNG figures only.

The transferred structural support is:

```text
L_REA_property = 204 pixels
```

With a property-map pixel size of \(0.01~\mathrm{mm}\), this corresponds to:

```text
2.04 mm
```

**Output directory:**

```text
property_full_homogenized_property_masks_STRICT/
```

The script writes:

```text
figure5_property_full_homogenized_property_masks_STRICT_kappa_app_curves.png
figure5_property_full_homogenized_property_masks_STRICT_C_app_curves.png
figure5_property_full_homogenized_property_masks_STRICT_E_app_curves.png
figure5_property_full_homogenized_property_masks_STRICT_summary.csv
figure5_property_full_homogenized_property_masks_STRICT_window_values.csv
```

The CSV files contain both the aggregated size-dependent results and the individual accepted-window calculations.

**Monitor the run:**

```bash
tail -f run_property_full_homogenized_property_masks_STRICT.log
```

---

### Workflow summary

```text
BSE TIFF images
      |
      v
detect_stationary_domains_general.py
      |
      v
BSE stationary masks
      |
      v
rev_stationary_mask_cli_v24_general.py
      |
      v
structural REA: 1536 BSE pixels, approximately 2.01 mm
      |
      | grid conversion
      v
predicted property-grid support: 204 pixels, approximately 2.04 mm
      |
      +--------------------------------------+
                                             |
QEMSCAN kappa, E, and nu CSV maps            |
      |                                      |
      v                                      |
detect_stationary_domains_general.py         |
      |                                      |
      v                                      |
property-grid stationary masks               |
      |                                      |
      +-------------------+------------------+
                          |
                          v
property_full_homogenized_inside_property_masks.py
                          |
                          v
conductivity, stiffness, and Young-modulus convergence curves
```

---

### Running long calculations with `nohup`

The BSE analysis and especially the property-level finite-element calculations can require substantial time and memory. The commands above use `nohup` so they continue after the terminal is closed.

The general form is:

```bash
nohup python script.py [arguments] > run.log 2>&1 &
```

Check the process:

```bash
ps -ef | grep python
```

Watch the log:

```bash
tail -f run.log
```

Stop a job by finding its process identifier and running:

```bash
kill PROCESS_ID
```

Use:

```bash
kill -9 PROCESS_ID
```

only when a normal `kill` does not stop the process.

---

### Important reproducibility notes

- Keep the ordering of the seven BSE images consistent throughout the structural analysis.
- Keep the ordering of `kappaMap_i.csv`, `EMap_i.csv`, and `nuMap_i.csv` identical.
- Property maps and their corresponding masks must have exactly the same dimensions.
- Use raw numerical property maps containing physical values.
- `--max-side 0` is essential for reproducing the full-resolution BSE analysis.
- The stationary masks are part of the calculation, not only visualization files.
- A candidate window is accepted only when its stationary-mask coverage satisfies `--min-mask-fraction`.
- The smoothed curves in the BSE summary figure are display aids. They do not change the numerical structural-support selection.
- `--no-pdf` avoids slow vector export of high-resolution image panels. Remove this option when PDF output is required.
- The property-level elasticity calculation is the most expensive stage. Reducing `--n-window-samples` or the largest tested window size is useful for preliminary tests, but manuscript reproduction should use the command reported above.
- The scripts use a deterministic random seed by default for property-window sampling unless `--rng-seed` is changed.

---

### Troubleshooting

#### A stationary mask cannot be found

Confirm that the image stem and mask suffix match. For example:

```text
Input image: 1_cut.tif
Expected mask: stationary_domains_FULLRES_R2048_STRICT/1_cut_stationary_mask.tif
```

#### Property mask naming does not match

Step 4 expects:

```text
stationary_property_combined_STRICT/property_combined_1_stationary_mask.tif
```

through:

```text
stationary_property_combined_STRICT/property_combined_7_stationary_mask.tif
```

#### Property arrays have different shapes

Each conductivity, Young's-modulus, Poisson-ratio, and mask array belonging to the same section must have identical dimensions.

#### No valid windows are accepted

Check:

- the mask itself;
- `--min-mask-fraction`;
- the requested window size;
- whether the stationary region is large enough to contain that support.

#### The elasticity calculation is slow

This is expected for large windows because the periodic Q4 finite-element system grows rapidly with the number of pixels. Test the workflow first with fewer support sizes or:

```text
--n-window-samples 1
```

For conductivity-only testing, the property script also supports:

```text
--skip-elastic
```

Do not use these reduced settings for the final manuscript reproduction unless they match the reported calculation.

#### Memory usage is too high

Large full-resolution images and large finite-element windows can require substantial memory. Run one stage at a time, monitor memory use, and avoid launching multiple property-homogenization jobs simultaneously.

---

### Associated manuscript

**Representativeness emerges from statistical admissibility and long-wavelength structural convergence**

The repository provides the analysis and homogenization scripts used for the two-dimensional numerical validation reported in the manuscript.

Citation information will be added after publication.

---

### Corresponding authors

Dr. Christian Tantardini  
Rice University  
`christian.tantardini@rice.edu`

Prof. Eduardo Garzanti  
University of Milano-Bicocca  
`eduardo.garzanti@unimib.it`
