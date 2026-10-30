# Compact figure inputs and derived tables

These files preserve data used in the revised figures without duplicating the
large BSE dataset. Original data and cell responses are from Zenodo version
10.5281/zenodo.21717613. The `.npz` files contain numerical arrays only and can
be opened with `numpy.load(path, allow_pickle=False)`.

| File | Contents and provenance |
|---|---|
| `figure2_plot_data.npz` | Eleven-support mean and spectral-residual series, SEMs, thresholds, selected/reference sizes, and selected ensemble radial spectra. Preserves the original numerical Figure 2 series. |
| `property_window_values.csv` | Copy of the archived 160 solved-window conductivity, stiffness, and axial-modulus responses. |
| `qemscan_image1_maps.npz` | The numerical conductivity, Young's modulus, and Poisson-ratio maps for image 1, losslessly repackaged from the archived CSVs. Used to reconstruct the supplementary map figure. |
| `phase_area_fractions.csv` | Counts and area fractions by image and phase, classified jointly from all three numerical property maps. |
| `phase_area_fractions_pooled.csv` | Area-weighted pooled phase counts and percentages across the seven complete images. |
| `axial_shear_window_values.csv` | `Ex`, `Ey`, and `Gxy` obtained from the inverse of each archived full 3-by-3 stiffness matrix. |
| `axial_shear_summary.csv` | Support-wise means and pooled-window SEMs for those moduli, with image/window counts. |
| `thermal_stiffness_summary.csv` | Support-wise means and pooled-window SEMs for archived conductivity and stiffness responses. |

Regenerate the CSV summaries and property figures with
`python figures/regenerate_revised_figures.py`. Supply `--data-root` to derive
phase counts and image-1 maps directly from the original extracted archive.
Figure 2 uses the same 2048-pixel reference and 1536-pixel selected support as
the archived baseline. Spectra shown in its plot are normalized for display;
this normalization does not enter the structural selection.

The 384-pixel property support has four windows from image 6 only. Supports
through 320 pixels have 12 windows from images 1, 5, and 6. Neither these
conditional domains nor the pooled-window error bars establish whole-sample
representativeness or a three-dimensional REV.
