# Compact figure inputs and derived tables

These files supply the numerical inputs for the manuscript figures without
duplicating the large BSE dataset. Image data and cell responses are from Zenodo version
10.5281/zenodo.21717613. The `.npz` files contain numerical arrays only and can
be opened with `numpy.load(path, allow_pickle=False)`.

| File | Contents and provenance |
|---|---|
| `figure2_plot_data.npz` | Eleven-support mean and spectral-residual series, SEMs, thresholds, selected/reference sizes, and selected ensemble radial spectra. Supplies the numerical Figure 2 series. |
| `property_window_values.csv` | Copy of the archived 160 solved-window conductivity, stiffness, and axial-modulus responses. |
| `qemscan_image1_maps.npz` | The numerical conductivity, Young's modulus, and Poisson-ratio maps for image 1, losslessly repackaged from the archived CSVs. Used to reconstruct the supplementary map figure. |
| `phase_area_fractions.csv` | Counts and area fractions by image and phase, classified jointly from all three numerical property maps. |
| `phase_area_fractions_pooled.csv` | Area-weighted pooled phase counts and percentages across the seven complete images. |
| `axial_shear_window_values.csv` | `Ex`, `Ey`, and `Gxy` obtained from the inverse of each archived full 3-by-3 stiffness matrix. |
| `axial_shear_summary.csv` | Support-wise means and pooled-window SEMs for those moduli, with image/window counts. |
| `thermal_stiffness_summary.csv` | Support-wise means and pooled-window SEMs for archived conductivity and stiffness responses. |

Generate the manuscript Figure 3 and CSV summaries from the bundled inputs with
`env -u REA_DATA_ROOT python figures/figure_3.py`. To use extracted study data,
supply `--data-root /path/to/data --study-preset`; the preset defines the study's
phase assignments, property-pixel calibration and support marker. Figure 2 is
assembled by `figures/figure_2.py` using `figure2_plot_data.npz` together with a
BSE image. See [figures/README.md](../figures/README.md) for the three manuscript
assembly commands and input paths.
Figure 2 uses the same 2048-pixel reference and 1536-pixel selected support as
the archived baseline. Spectra shown in its plot are normalized for display;
this normalization does not enter the structural selection.

The 384-pixel property support has four windows from image 6 only. Supports
through 320 pixels have 12 windows from images 1, 5, and 6. Neither these
conditional domains nor the pooled-window error bars establish whole-sample
representativeness or a three-dimensional REV.
