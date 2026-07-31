#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
REA convergence analysis for two-dimensional BSE gray-level images. In combined mode, multiple images are treated as independent fields of view and one publication-style summary figure is produced.

Physical interpretation
-----------------------
The BSE gray level is treated as an image-contrast field that reflects
mineralogical and textural heterogeneity. It is not interpreted as thermal
conductivity or as a direct quantitative composition map. Therefore the L-dependent scalar convergence panel
plots the centered L x L window average of the gray-intensity field:

    <Z_gray>_L = mean[Z_gray(x) over the L x L window].

The code also computes normalized convergence indicators for panel (c):
    1. eta_Z: five-point running-median plus monotone tail-maximum plateau envelope excluding the reference point of the ensemble-mean relative deviation of <Z_gray>_L from L_ref,
       computed sample by sample and then ensemble-averaged;
    2. eta_C: relative low-k second-order spectral residual with respect to L_ref;
    3. eta_ens: monotone tail envelope of the ensemble coefficient of variation across independent images.

No conductivity equation is solved in this version.

Numerical spectral definition
-----------------------------
The spectral descriptor is computed from the squared Fourier amplitude of the
mean-centred, Hann-tapered image field and then radially averaged. Thus, this
script evaluates the power spectrum directly. It does not independently
calculate a real-space autocovariance for a Wiener--Khinchin consistency test.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
Image.MAX_IMAGE_PIXELS = None

# Publication-style defaults: clean vector text, readable labels, consistent strokes.
plt.rcParams.update({
    "font.size": 9.0,
    "axes.labelsize": 9.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 7.0,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.minor.width": 0.6,
    "ytick.minor.width": 0.6,
    "lines.linewidth": 1.35,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# ============================================================
# Defaults
# ============================================================

ORIGINAL_PIXEL_SIZE = 1.0
MAX_SIDE_FOR_ANALYSIS = 1024
CROP_BOX = None
INVERT_GRAY = False
FIELD_MODE = "raw"          # raw gray level, or scaled to [Z_MIN, Z_MAX]
Z_MIN = 1.0
Z_MAX = 255.0
GRAY_LOW_PERCENTILE = 1.0
GRAY_HIGH_PERCENTILE = 99.0
WINDOW_SIZES = [64, 96, 128, 160, 192, 224, 256, 320, 384, 448, 512, 640, 768, 896]
REFERENCE_L = 896
TAU_Z = 2.5e-2
TAU_C = 3.0e-1
TAU_ENS = 1.0e-1
TAU_SHELL = 1.0e-1
SPECTRAL_BAND = "low"
LOW_K_FRACTION = 0.25
LOW_K_WEIGHT = "inverse"
ALPHA_SHELL = 0.80
N_COMMON_K_GRID = 300
TAIL_FRACTION = 0.50
LOG_FLOOR = 1.0e-4
SELECTED_L_FOR_SPECTRA = [64, 128, 256, 384, 512, 768, 896]
OUTPUT_PREFIX = "bse_gray_REA"
OUTPUT_DIR = Path(".")
SAVE_DPI = 500
N_WINDOW_SAMPLES = 49
ERRORBAR_MODE = "sem"
STATIONARY_MASK_DIR = None
STATIONARY_MASK_SUFFIX = "_stationary_mask.tif"
MIN_MASK_FRACTION = 0.98
# Panel-(b) display controls. The raw ensemble mean is always the source data;
# optional smoothing/reference band only improves visual readability and does
# not enter the REA decision.
PANEL_B_MODE = "raw+smooth"      # choices: raw, smooth, raw+smooth
PANEL_B_SMOOTH_WINDOW = 3          # odd running-median window in L-index space
PANEL_B_REFERENCE_BAND = True      # show +/- tau_Z band around reference mean
PANEL_C_MODE = "raw+smooth"      # choices: raw, smooth, raw+smooth
PANEL_C_SMOOTH_WINDOW = 3          # odd running-median window in L-index space


def parse_int_list(text: str) -> list[int]:
    if text is None or str(text).strip() == "":
        return []
    return [int(x.strip()) for x in str(text).split(",") if x.strip()]


def parse_crop_box(text: str | None):
    if text is None:
        return None
    if str(text).strip().lower() in {"none", "no", "false", "0"}:
        return None
    values = [int(x.strip()) for x in str(text).split(",") if x.strip()]
    if len(values) != 4:
        raise argparse.ArgumentTypeError(
            "Crop box must be 'left,upper,right,lower', for example '100,50,1800,1600'."
        )
    left, upper, right, lower = values
    if right <= left or lower <= upper:
        raise argparse.ArgumentTypeError("Invalid crop box: right must exceed left and lower must exceed upper.")
    return (left, upper, right, lower)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "REA convergence analysis for two-dimensional BSE gray-level images. "
            "The gray level is treated as an image-contrast field, not as conductivity."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("images", nargs="+", help="Input image file(s), e.g. image.tif image2.jpeg")
    parser.add_argument("--output-dir", default=".", help="Directory where figures are written.")
    parser.add_argument(
        "--output-prefix",
        default=None,
        help=(
            "Output prefix. For one image this is used exactly. For multiple images, "
            "the image stem is appended automatically. If omitted, each image stem is used."
        ),
    )

    parser.add_argument("--pixel-size", type=float, default=ORIGINAL_PIXEL_SIZE,
                        help="Pixel size of the original image in physical units. Use 1 for pixel units.")
    parser.add_argument("--max-side", type=int, default=MAX_SIDE_FOR_ANALYSIS,
                        help="Maximum side after downsampling. Use 0 to disable downsampling.")
    parser.add_argument("--crop", type=parse_crop_box, default=CROP_BOX,
                        help="Crop box in original image pixels: left,upper,right,lower. Use 'none' for no crop.")
    parser.add_argument("--invert-gray", action="store_true", default=INVERT_GRAY,
                        help="Use when brighter pixels should correspond to smaller Z_gray.")

    parser.add_argument("--field-mode", choices=["raw", "scaled"], default=FIELD_MODE,
                        help="Use raw gray intensity or robustly scale it to [z-min, z-max].")
    parser.add_argument("--z-min", type=float, default=Z_MIN,
                        help="Minimum value of the scaled gray-level analysis field.")
    parser.add_argument("--z-max", type=float, default=Z_MAX,
                        help="Maximum value of the scaled gray-level analysis field.")
    parser.add_argument("--gray-low-percentile", type=float, default=GRAY_LOW_PERCENTILE,
                        help="Lower robust percentile for scaled mapping.")
    parser.add_argument("--gray-high-percentile", type=float, default=GRAY_HIGH_PERCENTILE,
                        help="Upper robust percentile for scaled mapping.")

    parser.add_argument("--window-sizes", type=parse_int_list, default=WINDOW_SIZES,
                        help="Comma-separated window sizes, e.g. 64,96,128,256,512,896.")
    parser.add_argument("--reference-L", type=int, default=REFERENCE_L,
                        help="Requested reference window. If absent or invalid, the largest valid L is used.")
    parser.add_argument("--selected-spectra", type=parse_int_list, default=SELECTED_L_FOR_SPECTRA,
                        help="Comma-separated L values to show in the spectra figure.")

    parser.add_argument("--tau-Z", dest="tau_Z", type=float, default=TAU_Z,
                        help="Tolerance for relative change in L-window mean gray/Z value.")
    parser.add_argument("--tau-C", dest="tau_C", type=float, default=TAU_C,
                        help="Tolerance for the low-wavenumber second-order spectral residual.")
    parser.add_argument(
        "--tau-ens", "--tau-boundary", dest="tau_ens", type=float, default=TAU_ENS,
        help=(
            "Tolerance for the auxiliary ensemble-reproducibility residual. "
            "--tau-boundary is retained as a backward-compatible alias."
        ),
    )
    parser.add_argument(
        "--tau-shell", type=float, default=TAU_SHELL,
        help="Tolerance for the auxiliary boundary-shell residual in single-image diagnostic mode.",
    )
    parser.add_argument("--spectral-band", choices=["low", "full"], default=SPECTRAL_BAND,
                        help="Spectral metric for panel (c): low-k REA-sensitive band or full resolved band.")
    parser.add_argument("--low-k-fraction", type=float, default=LOW_K_FRACTION,
                        help="Fraction of the common k-grid used when --spectral-band low.")
    parser.add_argument("--low-k-weight", choices=["uniform", "inverse"], default=LOW_K_WEIGHT,
                        help="Weighting used for the low-k spectral residual.")
    parser.add_argument("--alpha-shell", type=float, default=ALPHA_SHELL,
                        help="Inner-window fraction used to define the outer boundary shell.")
    parser.add_argument("--n-common-k-grid", type=int, default=N_COMMON_K_GRID,
                        help="Number of points in the common spectral grid.")
    parser.add_argument("--tail-fraction", type=float, default=TAIL_FRACTION,
                        help="Persistent REA tail fraction.")
    parser.add_argument("--log-floor", type=float, default=LOG_FLOOR,
                        help="Numerical floor for semilog diagnostic plots.")
    parser.add_argument("--dpi", type=int, default=SAVE_DPI,
                        help="PNG resolution for the saved summary figure.")
    parser.add_argument("--n-window-samples", type=int, default=N_WINDOW_SAMPLES,
                        help="Approximate number of L x L windows sampled over the image for mean and error bars.")
    parser.add_argument("--errorbar-mode", choices=["std", "sem"], default=ERRORBAR_MODE,
                        help="Use standard deviation or standard error of the sampled L-window averages as error bars.")
    parser.add_argument("--save-diagnostics", action="store_true",
                        help="Also save the old separate diagnostic figures. By default only the 4-panel summary is saved.")
    parser.add_argument("--no-pdf", action="store_true",
                        help="Save PNG only. This is recommended for batch runs because PDF export of high-resolution image panels can be slow.")
    parser.add_argument("--combine-images", action="store_true", default=True,
                        help="Treat all input images as independent fields of view and produce one combined 4-panel figure. Default: enabled.")
    parser.add_argument("--separate-images", action="store_true",
                        help="Process each input image separately. Use only for diagnostics; manuscript mode is combined.")
    parser.add_argument("--representative-image", type=int, default=1,
                        help="1-based index of the input image shown in panel (a) and used for panel (d) spectra in combined mode.")
    parser.add_argument("--stationary-mask-dir", default=None,
                        help="Directory containing stationary masks named <image_stem><mask-suffix>. If omitted, full images are used.")
    parser.add_argument("--stationary-mask-suffix", default=STATIONARY_MASK_SUFFIX,
                        help="Suffix used to find stationary masks in --stationary-mask-dir.")
    parser.add_argument("--min-mask-fraction", type=float, default=MIN_MASK_FRACTION,
                        help="Minimum fraction of pixels inside stationary mask required for an L x L window to be used.")

    # Panel-(b) controls. These options change only the display of the ensemble
    # mean curve. They do not change the structural REA criterion, which remains
    # the low-k second-order spectral residual in panel (c).
    parser.add_argument(
        "--panel-b-mode",
        choices=["raw", "smooth", "raw+smooth"],
        default=PANEL_B_MODE,
        help=(
            "Display mode for panel (b). 'raw' plots the measured ensemble mean; "
            "'smooth' plots a running-median trend; 'raw+smooth' plots raw points "
            "faintly plus the running-median trend. This is for readability only."
        ),
    )
    parser.add_argument(
        "--panel-b-smooth-window",
        type=int,
        default=PANEL_B_SMOOTH_WINDOW,
        help=(
            "Odd running-median width, in number of tested L values, used only for "
            "the optional panel-(b) smoothed trend."
        ),
    )
    parser.add_argument(
        "--no-panel-b-reference-band",
        action="store_true",
        help=(
            "Do not draw the +/- tau_Z reference band around the finite reference "
            "mean in panel (b). The band is a visual guide only."
        ),
    )

    # Panel-(c) controls. These options modify only the visual presentation of
    # the structural residual curve. The numerical REA criterion is always
    # computed from the unsmoothed eta_C values written to the CSV/log.
    parser.add_argument(
        "--panel-c-mode",
        choices=["raw", "smooth", "raw+smooth"],
        default=PANEL_C_MODE,
        help=(
            "Display mode for panel (c). 'raw' plots the measured structural "
            "residual; 'smooth' plots a running-median trend; 'raw+smooth' "
            "shows faint raw points plus a smoothed trend. This is for visual "
            "readability only and does not change the REA criterion."
        ),
    )
    parser.add_argument(
        "--panel-c-smooth-window",
        type=int,
        default=PANEL_C_SMOOTH_WINDOW,
        help=(
            "Odd running-median width, in number of tested L values, used only "
            "for the optional panel-(c) smoothed structural residual."
        ),
    )
    return parser


def apply_args_to_globals(args):
    global ORIGINAL_PIXEL_SIZE, MAX_SIDE_FOR_ANALYSIS, CROP_BOX, INVERT_GRAY
    global FIELD_MODE, Z_MIN, Z_MAX, GRAY_LOW_PERCENTILE, GRAY_HIGH_PERCENTILE
    global WINDOW_SIZES, REFERENCE_L, TAU_Z, TAU_C, TAU_ENS, TAU_SHELL
    global SPECTRAL_BAND, LOW_K_FRACTION, LOW_K_WEIGHT
    global ALPHA_SHELL, N_COMMON_K_GRID, TAIL_FRACTION, LOG_FLOOR, SELECTED_L_FOR_SPECTRA
    global OUTPUT_DIR, SAVE_DPI, N_WINDOW_SAMPLES, ERRORBAR_MODE
    global STATIONARY_MASK_DIR, STATIONARY_MASK_SUFFIX, MIN_MASK_FRACTION
    global PANEL_B_MODE, PANEL_B_SMOOTH_WINDOW, PANEL_B_REFERENCE_BAND
    global PANEL_C_MODE, PANEL_C_SMOOTH_WINDOW

    ORIGINAL_PIXEL_SIZE = float(args.pixel_size)
    MAX_SIDE_FOR_ANALYSIS = None if args.max_side == 0 else int(args.max_side)
    CROP_BOX = args.crop
    INVERT_GRAY = bool(args.invert_gray)
    FIELD_MODE = args.field_mode
    Z_MIN = float(args.z_min)
    Z_MAX = float(args.z_max)
    GRAY_LOW_PERCENTILE = float(args.gray_low_percentile)
    GRAY_HIGH_PERCENTILE = float(args.gray_high_percentile)
    WINDOW_SIZES = [int(x) for x in args.window_sizes]
    REFERENCE_L = int(args.reference_L)
    TAU_Z = float(args.tau_Z)
    TAU_C = float(args.tau_C)
    TAU_ENS = float(args.tau_ens)
    TAU_SHELL = float(args.tau_shell)
    SPECTRAL_BAND = str(args.spectral_band)
    LOW_K_FRACTION = float(args.low_k_fraction)
    LOW_K_WEIGHT = str(args.low_k_weight)
    ALPHA_SHELL = float(args.alpha_shell)
    N_COMMON_K_GRID = int(args.n_common_k_grid)
    TAIL_FRACTION = float(args.tail_fraction)
    LOG_FLOOR = float(args.log_floor)
    SELECTED_L_FOR_SPECTRA = [int(x) for x in args.selected_spectra]
    OUTPUT_DIR = Path(args.output_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    SAVE_DPI = int(args.dpi)
    N_WINDOW_SAMPLES = int(args.n_window_samples)
    ERRORBAR_MODE = str(args.errorbar_mode)
    STATIONARY_MASK_DIR = Path(args.stationary_mask_dir) if args.stationary_mask_dir else None
    STATIONARY_MASK_SUFFIX = str(args.stationary_mask_suffix)
    MIN_MASK_FRACTION = float(args.min_mask_fraction)
    PANEL_B_MODE = str(args.panel_b_mode)
    PANEL_B_SMOOTH_WINDOW = int(args.panel_b_smooth_window)
    if PANEL_B_SMOOTH_WINDOW < 1:
        raise ValueError("--panel-b-smooth-window must be >= 1")
    if PANEL_B_SMOOTH_WINDOW % 2 == 0:
        PANEL_B_SMOOTH_WINDOW += 1
        print(f"Adjusted --panel-b-smooth-window to odd value: {PANEL_B_SMOOTH_WINDOW}")
    PANEL_B_REFERENCE_BAND = not bool(args.no_panel_b_reference_band)
    PANEL_C_MODE = str(args.panel_c_mode)
    PANEL_C_SMOOTH_WINDOW = int(args.panel_c_smooth_window)
    if PANEL_C_SMOOTH_WINDOW < 1:
        raise ValueError("--panel-c-smooth-window must be >= 1")
    if PANEL_C_SMOOTH_WINDOW % 2 == 0:
        PANEL_C_SMOOTH_WINDOW += 1
        print(f"Adjusted --panel-c-smooth-window to odd value: {PANEL_C_SMOOTH_WINDOW}")


# ============================================================
# Image and field construction
# ============================================================


def load_gray_image_as_array(path: str | os.PathLike, crop_box=None, max_side: int | None = 1024):
    Image.MAX_IMAGE_PIXELS = None
    img = Image.open(path)
    original_size = img.size

    if crop_box is not None:
        img = img.crop(crop_box)
    img = img.convert("L")
    pre_downsample_size = img.size

    scale = 1.0
    if max_side is not None and max(img.size) > max_side:
        scale = max_side / float(max(img.size))
        new_size = (
            max(1, int(round(img.size[0] * scale))),
            max(1, int(round(img.size[1] * scale))),
        )
        img = img.resize(new_size, Image.Resampling.LANCZOS)

    dx_eff = ORIGINAL_PIXEL_SIZE / scale
    gray = np.asarray(img, dtype=float)
    analysis_size = img.size

    print("\nLoaded image")
    print("=" * 72)
    print(f"Input path                    : {path}")
    print(f"Original size                 : {original_size[0]} x {original_size[1]} pixels")
    if crop_box is not None:
        print(f"Crop box                      : {crop_box}")
        print(f"Size after crop               : {pre_downsample_size[0]} x {pre_downsample_size[1]} pixels")
    print(f"Analysis size                 : {analysis_size[0]} x {analysis_size[1]} pixels")
    print(f"Downsample scale              : {scale:.6g}")
    print(f"Effective pixel size          : {dx_eff:.6g}")
    print(f"Gray range                    : {gray.min():.6g} to {gray.max():.6g}")
    print("=" * 72)
    return gray, dx_eff, original_size, analysis_size


def load_stationary_mask_for_image(image_path: str | os.PathLike, target_shape: tuple[int, int]) -> np.ndarray | None:
    """Load the stationary-domain mask corresponding to one image.

    The mask is expected at:
        <stationary-mask-dir>/<image_stem><stationary-mask-suffix>

    If the mask shape differs from the analysis image shape, it is resized with
    nearest-neighbour interpolation so that image and mask live on the same grid.
    """
    if STATIONARY_MASK_DIR is None:
        return None

    image_path = Path(image_path)
    mask_path = Path(STATIONARY_MASK_DIR) / f"{image_path.stem}{STATIONARY_MASK_SUFFIX}"
    if not mask_path.exists():
        raise FileNotFoundError(
            f"Stationary mask not found for image '{image_path.name}': {mask_path}\n"
            "Check --stationary-mask-dir and --stationary-mask-suffix."
        )

    m = Image.open(mask_path).convert("L")
    target_h, target_w = target_shape
    if m.size != (target_w, target_h):
        m = m.resize((target_w, target_h), Image.Resampling.NEAREST)

    mask = np.asarray(m, dtype=np.uint8) > 0
    frac = float(mask.mean())
    print(f"Loaded stationary mask        : {mask_path}")
    print(f"Mask shape                    : {mask.shape[1]} x {mask.shape[0]} pixels")
    print(f"Mask stationary pixel fraction: {frac:.6g}")
    return mask


def binary_integral_image(mask: np.ndarray) -> np.ndarray:
    """Integral image for a boolean stationary mask."""
    arr = mask.astype(np.int32)
    return np.pad(arr.cumsum(axis=0).cumsum(axis=1), ((1, 0), (1, 0)), mode="constant")


def integral_window_sum(ii: np.ndarray, y0: int, x0: int, y1: int, x1: int) -> int:
    return int(ii[y1, x1] - ii[y0, x1] - ii[y1, x0] + ii[y0, x0])


def masked_window_fraction(mask_ii: np.ndarray, y: int, x: int, L: int) -> float:
    return integral_window_sum(mask_ii, y, x, y + L, x + L) / float(L * L)


def valid_masked_starts(mask: np.ndarray | None, L: int, approx_n: int, min_fraction: float):
    """Return deterministic grid starts whose LxL window lies inside the stationary mask.

    If no mask is supplied, this returns the usual regular-grid starts over the
    full image. If a mask is supplied, starts are retained only if the fraction
    of stationary-mask pixels inside the LxL window is >= min_fraction.
    """
    if mask is None:
        H = W = None
        return None

    H, W = mask.shape
    if L > min(H, W):
        return []

    n_axis = max(1, int(np.ceil(np.sqrt(max(1, approx_n)))))
    ys = np.linspace(0, H - L, n_axis).round().astype(int)
    xs = np.linspace(0, W - L, n_axis).round().astype(int)
    starts = {(int(y), int(x)) for y in ys for x in xs}

    # Always try the centered window too.
    starts.add((H // 2 - L // 2, W // 2 - L // 2))

    ii = binary_integral_image(mask)
    out = []
    for y, x in sorted(starts):
        if masked_window_fraction(ii, y, x, L) >= min_fraction:
            out.append((y, x))

    return out


def find_representative_masked_window(arr: np.ndarray, L: int, mask: np.ndarray | None = None,
                                      approx_n: int = 49, min_fraction: float = 0.98):
    """Return one LxL field for spectral and boundary-shell diagnostics.

    The selected window is the valid stationary-mask window whose center is
    closest to the image center. If no mask is supplied, the centered window is
    returned. If no valid stationary window exists, None is returned.
    """
    H, W = arr.shape
    if L > min(H, W):
        return None

    if mask is None:
        return centered_square_window(arr, L)

    starts = valid_masked_starts(mask, L, approx_n=max(approx_n, N_WINDOW_SAMPLES), min_fraction=min_fraction)
    if not starts:
        return None

    cy0 = H / 2.0
    cx0 = W / 2.0
    y, x = min(starts, key=lambda yx: ((yx[0] + 0.5 * L - cy0) ** 2 + (yx[1] + 0.5 * L - cx0) ** 2))
    return arr[y:y + L, x:x + L].copy()


def has_valid_stationary_window(mask: np.ndarray | None, L: int, min_fraction: float) -> bool:
    if mask is None:
        return True
    starts = valid_masked_starts(mask, L, approx_n=49, min_fraction=min_fraction)
    return len(starts) > 0


def gray_to_analysis_field(gray: np.ndarray) -> np.ndarray:
    """Convert a grayscale image to the analysis field Z_gray(x)."""
    g = gray.astype(float)
    if INVERT_GRAY:
        g = 255.0 - g

    if FIELD_MODE.lower() == "raw":
        return g

    if FIELD_MODE.lower() != "scaled":
        raise ValueError("FIELD_MODE must be either 'raw' or 'scaled'.")

    g_low = np.percentile(g, GRAY_LOW_PERCENTILE)
    g_high = np.percentile(g, GRAY_HIGH_PERCENTILE)
    if not np.isfinite(g_low) or not np.isfinite(g_high) or g_high <= g_low:
        raise ValueError("Invalid gray percentile range for scaling.")

    gn = (g - g_low) / (g_high - g_low)
    gn = np.clip(gn, 0.0, 1.0)
    return Z_MIN + gn * (Z_MAX - Z_MIN)


def centered_square_window(arr: np.ndarray, L: int) -> np.ndarray:
    H, W = arr.shape
    if L > min(H, W):
        raise ValueError(f"Window L={L} exceeds min image dimension {min(H, W)}.")
    i0 = H // 2 - L // 2
    j0 = W // 2 - L // 2
    return arr[i0:i0 + L, j0:j0 + L].copy()




def centered_window_bounds(shape: tuple[int, int], L: int) -> tuple[int, int, int, int]:
    H, W = shape
    i0 = H // 2 - L // 2
    j0 = W // 2 - L // 2
    return j0, i0, L, L


def sample_square_windows(arr: np.ndarray, L: int, approx_n: int = 49,
                          mask: np.ndarray | None = None,
                          min_mask_fraction: float = 0.98) -> list[np.ndarray]:
    """Sample L x L windows on a regular grid.

    If a stationary-domain mask is provided, only windows whose fraction inside
    the mask is at least min_mask_fraction are retained. This makes panel (b)
    and the ensemble statistics conditional on the stationary domain.
    """
    H, W = arr.shape
    if L > min(H, W):
        raise ValueError(f"Window L={L} exceeds min image dimension {min(H, W)}.")

    if mask is None:
        if approx_n <= 1:
            return [centered_square_window(arr, L)]

        n_axis = max(1, int(np.ceil(np.sqrt(approx_n))))
        ys = np.linspace(0, H - L, n_axis).round().astype(int)
        xs = np.linspace(0, W - L, n_axis).round().astype(int)
        starts = {(int(y), int(x)) for y in ys for x in xs}
        starts.add((H // 2 - L // 2, W // 2 - L // 2))
    else:
        starts = set(valid_masked_starts(mask, L, approx_n=approx_n, min_fraction=min_mask_fraction))

    windows = [arr[y:y + L, x:x + L].copy() for y, x in sorted(starts)]
    return windows


def window_mean_statistics(arr: np.ndarray, L: int, approx_n: int = 49,
                           mask: np.ndarray | None = None,
                           min_mask_fraction: float = 0.98):
    windows = sample_square_windows(arr, L, approx_n, mask=mask, min_mask_fraction=min_mask_fraction)
    if len(windows) == 0:
        return np.nan, np.nan, np.nan, 0
    vals = np.array([float(np.mean(w)) for w in windows], dtype=float)
    mean = float(np.mean(vals))
    std = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0
    sem = float(std / np.sqrt(len(vals))) if len(vals) > 1 else 0.0
    return mean, std, sem, int(len(vals))

def validate_window_sizes(image_shape: tuple[int, int]) -> list[int]:
    H, W = image_shape
    max_L = min(H, W)
    valid = sorted({int(L) for L in WINDOW_SIZES if int(L) <= max_L})
    if len(valid) < 4:
        raise ValueError(
            f"Too few valid window sizes for image shape {image_shape}. "
            f"Largest possible L is {max_L}. Adjust --window-sizes or --max-side."
        )
    if REFERENCE_L not in valid:
        print(
            f"Warning: requested reference-L={REFERENCE_L} is not valid for this image. "
            f"Using largest valid L={valid[-1]} as reference."
        )
    return valid


# ============================================================
# Spectral and auxiliary boundary-shell diagnostics
# ============================================================


def radial_average_spectrum(field: np.ndarray, dx: float = 1.0):
    """Radially averaged FFT power spectrum of the mean-centred Z_gray field.

    This version uses vectorized binning with np.bincount. It is numerically
    equivalent to the original loop but much faster for batch processing.
    """
    f = field.astype(float) - np.mean(field)
    L = field.shape[0]
    taper_1d = np.hanning(L)
    f = f * np.outer(taper_1d, taper_1d)

    F = np.fft.fftshift(np.fft.fft2(f))
    P = np.abs(F) ** 2 / field.size

    ky = np.fft.fftshift(np.fft.fftfreq(L, d=dx)) * 2.0 * np.pi
    kx = np.fft.fftshift(np.fft.fftfreq(L, d=dx)) * 2.0 * np.pi
    KX, KY = np.meshgrid(kx, ky, indexing="xy")
    KR = np.sqrt(KX * KX + KY * KY)

    kr = KR.ravel()
    pr = P.ravel()
    nbins = max(16, L // 4)
    bins = np.linspace(0.0, kr.max(), nbins + 1)
    which = np.searchsorted(bins, kr, side="right") - 1
    valid = (which >= 0) & (which < nbins)

    sums = np.bincount(which[valid], weights=pr[valid], minlength=nbins)
    counts = np.bincount(which[valid], minlength=nbins).astype(float)
    mask = counts > 0

    k_centers = 0.5 * (bins[:-1] + bins[1:])
    p_radial = np.zeros(nbins, dtype=float)
    p_radial[mask] = sums[mask] / counts[mask]
    return k_centers[mask], p_radial[mask]


def trapezoid_weights(x: np.ndarray):
    x = np.asarray(x)
    w = np.zeros_like(x)
    if len(x) == 1:
        w[0] = 1.0
        return w
    w[0] = 0.5 * (x[1] - x[0])
    w[-1] = 0.5 * (x[-1] - x[-2])
    if len(x) > 2:
        w[1:-1] = 0.5 * (x[2:] - x[:-2])
    return w


def interpolate_spectrum(k, p, q_grid):
    return np.interp(q_grid, k, p)


def low_k_spectral_residual(k, p, k_ref, p_ref, q_grid, weights):
    """Relative spectral residual used for panel (c).

    The default is a low-k metric because REA convergence is controlled by
    large-scale heterogeneity. The residual is always normalized by the
    reference spectrum, so it is dimensionless and comparable across samples:

        eta_C(L) = ||C_L - C_ref|| / (||C_ref|| + eps).

    With --spectral-band low, only k <= LOW_K_FRACTION * kmax is used. With
    --low-k-weight inverse, small wavenumbers receive slightly larger weight.
    """
    C = interpolate_spectrum(k, p, q_grid)
    C_ref = interpolate_spectrum(k_ref, p_ref, q_grid)

    valid = np.isfinite(C) & np.isfinite(C_ref)
    if SPECTRAL_BAND == "low":
        kmax = float(np.nanmax(q_grid))
        kc = max(float(LOW_K_FRACTION) * kmax, np.nanmin(q_grid[q_grid > 0]) if np.any(q_grid > 0) else 0.0)
        valid &= q_grid <= kc

    if np.count_nonzero(valid) < 3:
        return np.nan

    w = np.array(weights, dtype=float)
    if SPECTRAL_BAND == "low" and LOW_K_WEIGHT == "inverse":
        positive = q_grid[q_grid > 0]
        eps_k = float(np.nanmin(positive)) if positive.size else 1.0
        w = w / (q_grid + eps_k)

    numerator = np.sum(w[valid] * (C[valid] - C_ref[valid]) ** 2)
    denominator = np.sum(w[valid] * C_ref[valid] ** 2) + 1.0e-30
    return float(np.sqrt(numerator / denominator))


# Backward-compatible name used by older branches in this file.
def low_k_spectral_residual_full_band(k, p, k_ref, p_ref, q_grid, weights):
    return low_k_spectral_residual(k, p, k_ref, p_ref, q_grid, weights)


def robust_material_scale(field: np.ndarray, p_low: float = 5.0, p_high: float = 95.0) -> float:
    """Robust gray-level contrast scale used only by the auxiliary boundary-shell metric."""
    finite_vals = field[np.isfinite(field)]
    if finite_vals.size == 0:
        return np.nan
    scale = float(np.percentile(finite_vals, p_high) - np.percentile(finite_vals, p_low))
    return scale


def boundary_shell_diagnostic(field: np.ndarray, alpha: float = 0.80, material_scale: float | None = None):
    """Compare outer-shell mean of Z_gray with full-window mean.

    Auxiliary single-image boundary-shell consistency residual,
    eta = |<Z>_shell - <Z>_full| / (S_Z + eps),
    where S_Z is a fixed robust material-contrast scale for the whole image/sample.

    This changes only eta_boundary.  mean_full and mean_shell are returned unchanged.
    """
    L = field.shape[0]
    inner = int(np.floor(alpha * L))
    inner = max(2, min(inner, L - 2))
    start = L // 2 - inner // 2
    end = start + inner

    shell_mask = np.ones_like(field, dtype=bool)
    shell_mask[start:end, start:end] = False

    mean_full = float(np.mean(field))
    mean_shell = float(np.mean(field[shell_mask]))
    eps = 1.0e-30

    if material_scale is None:
        material_scale = robust_material_scale(field)

    if not np.isfinite(material_scale) or material_scale <= 0.0:
        eta = np.nan
    else:
        eta = abs(mean_shell - mean_full) / (material_scale + eps)
    return eta, mean_full, mean_shell


# ============================================================
# REA selection
# ============================================================


def check_window_ok(r, tau_Z, tau_C):
    if np.isnan(r["delta_Z"]):
        return False
    if r["is_reference"]:
        return False
    if np.isnan(r["eta_C"]):
        return False
    return (
        r["delta_Z"] <= tau_Z
        and r["eta_C"] <= tau_C
    )


def find_persistent_rea(results, tau_Z, tau_C, tail_fraction=0.50):
    """Return the first non-reference L whose larger non-reference tail remains converged.

    The primary structural REA criterion uses only the persistent apparent-mean
    residual and the low-wavenumber spectral residual. Auxiliary ensemble and
    boundary-shell diagnostics are reported separately and do not enter this
    selection rule.
    """
    for i, r in enumerate(results):
        if np.isnan(r["delta_Z"]):
            continue
        if r["is_reference"]:
            continue
        tail = results[i:]
        tail_test = [rr for rr in tail if not rr["is_reference"]]
        if len(tail_test) == 0:
            continue
        required_tail_length = max(2, int(np.ceil(tail_fraction * len(tail))))
        if len(tail_test) < required_tail_length:
            continue
        if all(check_window_ok(rr, tau_Z, tau_C) for rr in tail_test):
            return r["L"]
    return None


# ============================================================
# Plot helpers
# ============================================================



def running_nanmedian_1d(values: np.ndarray, width: int) -> np.ndarray:
    """Return a local running median while ignoring NaNs.

    This is used only for visual smoothing of manuscript panels. It never
    enters the structural REA decision. The width is counted in tested-window
    indices, not in pixels.
    """
    values = np.asarray(values, dtype=float)
    width = max(1, int(width))
    if width % 2 == 0:
        width += 1
    half = width // 2
    out = np.full_like(values, np.nan, dtype=float)
    for j in range(values.size):
        i0 = max(0, j - half)
        i1 = min(values.size, j + half + 1)
        w = values[i0:i1]
        w = w[np.isfinite(w)]
        if w.size:
            out[j] = float(np.nanmedian(w))
    return out


def running_log_median_positive(values: np.ndarray, width: int) -> np.ndarray:
    """Running median in log10 space for positive residual curves.

    Panel (c) is plotted on a logarithmic y-axis. Smoothing directly in linear
    space can overweight isolated large residuals, while smoothing in log space
    preserves multiplicative structure. This function is display-only.
    """
    values = np.asarray(values, dtype=float)
    out = np.full_like(values, np.nan, dtype=float)
    pos = np.isfinite(values) & (values > 0.0)
    if not np.any(pos):
        return out
    logs = np.full_like(values, np.nan, dtype=float)
    logs[pos] = np.log10(np.maximum(values[pos], LOG_FLOOR))
    sm = running_nanmedian_1d(logs, width)
    finite = np.isfinite(sm)
    out[finite] = 10.0 ** sm[finite]
    return out


def diagnostic_array(results, key, hide_reference=False):
    values = []
    for r in results:
        if hide_reference and r["is_reference"]:
            values.append(np.nan)
            continue
        val = r[key]
        if np.isnan(val):
            values.append(np.nan)
        else:
            values.append(max(val, LOG_FLOOR))
    return np.array(values)


def add_rea_line(ax, accepted_L):
    if accepted_L is not None:
        ax.axvline(
            accepted_L,
            color="black",
            linestyle="--",
            linewidth=1.5,
            alpha=0.90,
            label=rf"$L_{{\rm REA}}={accepted_L}$",
        )


def save_png_pdf(fig, basename, dpi=None, save_pdf=True):
    if dpi is None:
        dpi = SAVE_DPI
    fig.savefig(f"{basename}.png", dpi=dpi, bbox_inches="tight")
    if save_pdf:
        fig.savefig(f"{basename}.pdf", bbox_inches="tight")


# ============================================================
# Plotting
# ============================================================


def make_figures(z_full, gray_full, results, accepted_L, kmax, save_diagnostics=False, no_pdf=False):
    Lvals = np.array([r["L"] for r in results])
    center_mean_Z = np.array([r["mean_Z"] for r in results])
    sample_mean_Z = np.array([r["sample_mean_Z"] for r in results])
    sample_std_Z = np.array([r["sample_std_Z"] for r in results])
    sample_sem_Z = np.array([r["sample_sem_Z"] for r in results])
    sample_n = np.array([r["sample_n"] for r in results])
    yerr = sample_sem_Z if ERRORBAR_MODE == "sem" else sample_std_Z

    raw_delta = diagnostic_array(results, "delta_Z", hide_reference=True)
    raw_etaC = diagnostic_array(results, "eta_C", hide_reference=True)
    raw_etab = diagnostic_array(results, "eta_boundary", hide_reference=False)

    crit_delta = raw_delta / TAU_Z
    crit_etaC = raw_etaC / TAU_C
    crit_etab = raw_etab / TAU_SHELL

    # The single-image workflow has no ensemble uncertainty estimate. Keep
    # zero-length error bars so this diagnostic mode remains executable.
    crit_etaC_err = np.zeros_like(crit_etaC, dtype=float)
    crit_etab_err = np.zeros_like(crit_etab, dtype=float)

    ymax = max(np.nanmax(crit_delta), np.nanmax(crit_etaC), np.nanmax(crit_etab), 1.5)
    ymin = LOG_FLOOR

    # Single manuscript-style 4-panel figure. Titles are intentionally omitted;
    # only panel labels (a)--(d) are kept, following the coauthor request.
    fig, axes = plt.subplots(2, 2, figsize=(7.85, 6.35))

    ax = axes[0, 0]
    im = ax.imshow(z_full, cmap="gray", origin="lower")
    # Draw accepted REA square if available; otherwise draw the reference square.
    box_L = accepted_L if accepted_L is not None else int(Lvals[-1])
    x0, y0, w, h = centered_window_bounds(z_full.shape, int(box_L))
    rect = plt.Rectangle((x0, y0), w, h, fill=False, edgecolor="white", linewidth=1.5)
    ax.add_patch(rect)
    ax.text(0.03, 0.95, "(a)", transform=ax.transAxes, ha="left", va="top",
            fontsize=11, fontweight="bold", color="white",
            bbox=dict(facecolor="black", alpha=0.35, edgecolor="none", pad=2.0))
    ax.set_xlabel("x [pixels]")
    ax.set_ylabel("y [pixels]")
    # Keep the gray-level colorbar, but place it with enough padding so it
    # does not overlap the x-axis tick labels in the compact 4-panel layout.
    cbar = fig.colorbar(
        im,
        ax=ax,
        orientation="vertical",
        fraction=0.035,
        pad=0.055,
        shrink=0.78,
    )
    cbar.set_label(r"$Z_{\rm gray}$", fontsize=8)
    cbar.ax.tick_params(labelsize=7, pad=1.5)

    ax = axes[0, 1]
    ax.plot(Lvals, sample_mean_Z, "o-", linewidth=1.4,
            markersize=3.8, label=r"sample mean")
    ax.plot(Lvals, center_mean_Z, "s--", linewidth=1.0, markersize=3.0,
            label=r"centered window")
    add_rea_line(ax, accepted_L)
    ax.text(0.03, 0.95, "(b)", transform=ax.transAxes, ha="left", va="top",
            fontsize=11, fontweight="bold")
    ax.set_xlabel(r"window size $L$ [pixels]")
    ax.set_ylabel(r"$\langle Z_{\rm gray}\rangle_L$")
    ax.grid(True, alpha=0.22)
    ax.legend(frameon=False, fontsize=7.5)

    ax = axes[1, 0]
    panel_c_ymin = max(1.0e-2, ymin)
    ax.fill_between([Lvals.min(), Lvals.max()], panel_c_ymin, 1.0, color="0.92", alpha=0.70, zorder=0)

    # Plot the full convergence curves, but show error bars only from the
    # accepted REA size onward. This avoids oversized pre-convergence
    # uncertainty bars while retaining the mean trend before L_REA.
    if accepted_L is None:
        err_mask = np.ones_like(Lvals, dtype=bool)
    else:
        err_mask = Lvals >= float(accepted_L)

    ax.plot(Lvals, crit_delta, "o-", linewidth=1.25, markersize=3.5,
            label=r"$\eta_Z^{\rm med\text{-}tail}/\tau_Z$")
    ax.plot(Lvals, crit_etaC, "s-", linewidth=1.25, markersize=3.5,
            label=r"$\eta_{\widehat C}/\tau_{\widehat C}$")
    ax.plot(Lvals, crit_etab, "^-", linewidth=1.25, markersize=3.5,
            label=r"$\eta_{\rm shell}/\tau_{\rm shell}$")

    # Blue is the five-point running-median plus monotone tail-maximum
    # convergence plateau and is plotted without an error bar. Orange keeps
    # its uncertainty bar. Green is the monotone tail envelope of the
    # auxiliary boundary-shell residual. The zero-length error bars
    # simply preserve a common plotting path with the ensemble workflow.
    ax.errorbar(Lvals[err_mask], crit_etaC[err_mask], yerr=crit_etaC_err[err_mask],
                fmt="none", ecolor="C1", capsize=2.0, linewidth=1.0)
    ax.errorbar(Lvals[err_mask], crit_etab[err_mask], yerr=crit_etab_err[err_mask],
                fmt="none", ecolor="C2", capsize=2.0, linewidth=1.0)
    ax.set_yscale("log")
    ax.axhline(1.0, color="black", linestyle=":", linewidth=1.1)
    add_rea_line(ax, accepted_L)
    ax.text(0.03, 0.95, "(c)", transform=ax.transAxes, ha="left", va="top",
            fontsize=11, fontweight="bold")
    ax.set_ylim(panel_c_ymin, ymax * 1.12)
    ax.set_xlabel(r"window size $L$ [pixels]")
    ax.set_ylabel("residual / tolerance")
    # Ensure that all available post-REA points remain visible, with a small
    # right margin so the last post-REA/reference point is not clipped by the
    # axis frame.
    if np.any(np.isfinite(Lvals)):
        xmin = float(np.nanmin(Lvals))
        xmax = float(np.nanmax(Lvals))
        pad = 0.06 * max(xmax - xmin, 1.0)
        ax.set_xlim(xmin - 0.02 * max(xmax - xmin, 1.0), xmax + pad)
    ax.grid(True, which="major", alpha=0.24)
    ax.grid(True, which="minor", alpha=0.10)
    ax.legend(frameon=False, fontsize=7.0, handlelength=2.0)

    ax = axes[1, 1]
    d_xvals = []
    d_yvals = []
    shown = set(SELECTED_L_FOR_SPECTRA)
    for r in results:
        if r["L"] in shown:
            k = r["spectrum_k"]
            p = r["spectrum_p"]
            mask = (k > 0) & (p > 0)
            if np.any(mask):
                y = p[mask] / (np.nanmax(p[mask]) + 1.0e-30)
                d_xvals.append(k[mask])
                d_yvals.append(y)
                ax.loglog(k[mask], y, linewidth=1.35, label=rf"$L={r['L']}$")
    # The k_REA marker was intentionally removed from panel (d) to avoid
    # cluttering the low-k spectral-collapse plot and its legend.

    # Dynamic log-axis limits so the selected spectra are fully visible.
    if d_xvals and d_yvals:
        dx = np.concatenate([np.asarray(v, dtype=float) for v in d_xvals])
        dy = np.concatenate([np.asarray(v, dtype=float) for v in d_yvals])
        dx = dx[np.isfinite(dx) & (dx > 0)]
        dy = dy[np.isfinite(dy) & (dy > 0)]
        if dx.size:
            ax.set_xlim(dx.min() / 1.15, dx.max() * 1.15)
        if dy.size:
            ax.set_ylim(max(dy.min() / 1.8, 1.0e-8), min(dy.max() * 1.45, 2.0))
    ax.text(0.03, 0.95, "(d)", transform=ax.transAxes, ha="left", va="top",
            fontsize=11, fontweight="bold")
    ax.set_xlabel(r"$\log_{10}(k\,[\mathrm{pixel}^{-1}])$")
    ax.set_ylabel(r"$\log_{10}\!\left(\widehat{C}_L(k)\right)$")
    ax.grid(True, which="both", alpha=0.22)
    ax.legend(frameon=False, fontsize=6.6, ncol=2, handlelength=1.8, columnspacing=0.9)

    fig.tight_layout(pad=0.9, w_pad=1.0, h_pad=0.9)
    save_png_pdf(fig, f"{OUTPUT_PREFIX}_4panels_style_summary", dpi=SAVE_DPI, save_pdf=not no_pdf)

    if save_diagnostics:
        # Optional legacy diagnostics, disabled by default to avoid many plots.
        fig2, ax = plt.subplots(figsize=(6.8, 4.8))
        ax.errorbar(Lvals, sample_mean_Z, yerr=yerr, fmt="o-", linewidth=2.0,
                    capsize=3, label=rf"mean ± {ERRORBAR_MODE}")
        ax.plot(Lvals, center_mean_Z, "s--", linewidth=1.5, label="centered window")
        add_rea_line(ax, accepted_L)
        ax.set_xlabel(r"window size $L$ [pixels]")
        ax.set_ylabel(r"$\langle Z_{\rm gray}\rangle_L$")
        ax.grid(True, alpha=0.25)
        ax.legend(frameon=False)
        save_png_pdf(fig2, f"{OUTPUT_PREFIX}_Zgray_mean_convergence", save_pdf=not no_pdf)

    plt.close("all")

    print("\nSaved figures:")
    if no_pdf:
        print(f"  {OUTPUT_PREFIX}_4panels_style_summary.png")
    else:
        print(f"  {OUTPUT_PREFIX}_4panels_style_summary.png / .pdf")
    if save_diagnostics:
        if no_pdf:
            print(f"  {OUTPUT_PREFIX}_Zgray_mean_convergence.png")
        else:
            print(f"  {OUTPUT_PREFIX}_Zgray_mean_convergence.png / .pdf")




# ============================================================
# Combined multi-image / ensemble mode
# ============================================================


def common_valid_window_sizes(shapes: list[tuple[int, int]]) -> list[int]:
    """Return window sizes valid for every image in the ensemble."""
    max_common_L = min(min(s) for s in shapes)
    valid = sorted({int(L) for L in WINDOW_SIZES if int(L) <= max_common_L})
    if len(valid) < 4:
        raise ValueError(
            f"Too few common valid window sizes for image shapes {shapes}. "
            f"Largest common possible L is {max_common_L}. Adjust --window-sizes or --max-side."
        )
    return valid


def common_valid_window_sizes_with_masks(shapes: list[tuple[int, int]], masks: list[np.ndarray | None]) -> list[int]:
    """Return window sizes valid for all images and all stationary masks.

    A window size is retained only if every image can provide at least one
    LxL window satisfying the stationary-mask fraction criterion.
    """
    valid = common_valid_window_sizes(shapes)
    if all(m is None for m in masks):
        return valid

    filtered = []
    for L in valid:
        ok_all = True
        for m in masks:
            if not has_valid_stationary_window(m, L, MIN_MASK_FRACTION):
                ok_all = False
                break
        if ok_all:
            filtered.append(L)

    if len(filtered) < 4:
        raise ValueError(
            "Too few window sizes remain after applying stationary-domain masks. "
            "Use smaller --window-sizes, decrease --min-mask-fraction, or generate larger masks."
        )
    if REFERENCE_L not in filtered:
        print(
            f"Warning: requested reference-L={REFERENCE_L} is not valid inside all stationary masks. "
            f"Using largest valid L={filtered[-1]} as reference."
        )
    return filtered


def compute_image_results(z_full: np.ndarray, dx_eff: float, valid_window_sizes: list[int], reference_L: int, stationary_mask: np.ndarray | None = None):
    """Compute all per-image quantities without plotting.

    Panel-(c) scalar convergence is now referenced to the same image's
    reference window L_ref, not to the previous L. This reduces artificial
    image-to-image spread and gives a dimensionless metric:

        eta_Z_i(L) = |<Z>_{L,i} - <Z>_{L_ref,i}| / (|<Z>_{L_ref,i}| + eps).

    The spectral metric eta_C is also relative to the same image's L_ref
    spectrum, optionally restricted to the low-k band.
    """
    kmax = np.pi / dx_eff
    q_grid = np.linspace(0.0, kmax, N_COMMON_K_GRID)
    weights = trapezoid_weights(q_grid)

    # Fixed sample-level scale used only for the green boundary-shell metric.
    # This prevents the green residual from being dominated by an L-dependent
    # denominator while preserving the boundary-vs-full-window numerator.
    boundary_material_scale = robust_material_scale(z_full)

    spectra = {}
    window_stats = {}
    boundary_stats = {}

    for L in valid_window_sizes:
        z_L = find_representative_masked_window(
            z_full, L, mask=stationary_mask,
            approx_n=N_WINDOW_SAMPLES, min_fraction=MIN_MASK_FRACTION
        )
        if z_L is None:
            spectra[L] = (np.array([np.nan]), np.array([np.nan]))
            boundary_stats[L] = (np.nan, np.nan, np.nan)
        else:
            spectra[L] = radial_average_spectrum(z_L, dx=dx_eff)
            boundary_stats[L] = boundary_shell_diagnostic(
                z_L, alpha=ALPHA_SHELL, material_scale=boundary_material_scale
            )
        window_stats[L] = window_mean_statistics(
            z_full, L, approx_n=N_WINDOW_SAMPLES,
            mask=stationary_mask, min_mask_fraction=MIN_MASK_FRACTION
        )

    k_ref, p_ref = spectra[reference_L]
    ref_sample_mean_Z = window_stats[reference_L][0]
    ref_center_mean_Z = boundary_stats[reference_L][1]

    results = []
    for L in valid_window_sizes:
        k, p = spectra[L]
        eta_bnd, center_mean_Z, mean_shell = boundary_stats[L]
        sample_mean_Z, sample_std_Z, sample_sem_Z, sample_n = window_stats[L]
        is_reference = L == reference_L

        eta_C = np.nan if is_reference else low_k_spectral_residual(
            k, p, k_ref, p_ref, q_grid, weights
        )

        # Use the sampled L-window mean for the composition convergence metric,
        # because panel (b) is also based on sampled L-window averages.
        eta_Z = abs(sample_mean_Z - ref_sample_mean_Z) / max(abs(ref_sample_mean_Z), 1.0e-30)
        # Extra diagnostic retained in the CSV/objects, based on the centered window.
        eta_Z_centered = abs(center_mean_Z - ref_center_mean_Z) / max(abs(ref_center_mean_Z), 1.0e-30)

        results.append({
            "L": L,
            "mean_Z": center_mean_Z,
            "mean_shell": mean_shell,
            "sample_mean_Z": sample_mean_Z,
            "sample_std_Z": sample_std_Z,
            "sample_sem_Z": sample_sem_Z,
            "sample_n": sample_n,
            "delta_Z": eta_Z,
            "eta_Z": eta_Z,
            "eta_Z_centered": eta_Z_centered,
            "eta_C": eta_C,
            "eta_boundary": eta_bnd,
            "spectrum_k": k,
            "spectrum_p": p,
            "is_reference": is_reference,
        })
    return results, kmax


def aggregate_ensemble_results(per_image_results: list[list[dict]], valid_window_sizes: list[int], reference_L: int):
    """Aggregate 7 images as independent fields of view and return one results list.

    Panel (b) uses the ensemble mean of the sampled window means.
    Error bars are image-to-image standard deviation or SEM, depending on --errorbar-mode.
    The combined analysis computes three dimensionless diagnostics. The persistent
    mean residual and low-wavenumber spectral residual define the structural
    REA. The ensemble coefficient of variation is retained as an auxiliary
    reproducibility check. The blue diagnostic is the systematic convergence
    of the ensemble mean gray-level statistic
    at size L relative to the ensemble mean at L_ref. Its error bar
    is a jackknife uncertainty of that scalar residual, not the raw spread of
    the seven images.
    Green is the seven-image image-to-image coefficient of variation of the L-window mean and is used
    as the auxiliary REA reproducibility residual.
    Panel (d) is generated from representative/ensemble-average spectra elsewhere.
    """
    n_images = len(per_image_results)
    agg = []

    # BLUE metric in panel (c): ensemble mean gray-level convergence/bias.
    # Use the ensemble mean at L_ref as the reference, so the blue residual
    # measures systematic finite-size bias of the apparent mean gray level.
    # The raw image-to-image spread is already represented by the green
    # ensemble reproducibility metric, so the blue curve is plotted without
    # an error bar in panel (c).
    ref_idx = valid_window_sizes.index(reference_L)
    ref_rows = [res[ref_idx] for res in per_image_results]
    ref_sample_vals = np.array([r["sample_mean_Z"] for r in ref_rows], dtype=float)
    ref_ensemble_mean_Z = float(np.nanmean(ref_sample_vals))

    for idx, L in enumerate(valid_window_sizes):
        rows = [res[idx] for res in per_image_results]
        center_vals = np.array([r["mean_Z"] for r in rows], dtype=float)
        sample_vals = np.array([r["sample_mean_Z"] for r in rows], dtype=float)
        eta_b_vals = np.array([r["eta_boundary"] for r in rows], dtype=float)
        eta_C_vals = np.array([r["eta_C"] for r in rows], dtype=float)
        delta_vals = np.array([r["delta_Z"] for r in rows], dtype=float)

        center_mean = float(np.nanmean(center_vals))
        sample_mean = float(np.nanmean(sample_vals))
        sample_std = float(np.nanstd(sample_vals, ddof=1)) if n_images > 1 else 0.0
        sample_sem = float(sample_std / np.sqrt(n_images)) if n_images > 1 else 0.0
        # GREEN metric in panel (c), ensemble REA reproducibility criterion.
        # The image-to-image spread among the independent fields of view is the diagnostic itself:
        #
        #     eta_ens(L) = std_i[ <Z>_L^{(i)} ] / ( |mean_i[ <Z>_L^{(i)} ]| + eps )
        #
        # The green error bar is not the original image-to-image spread again.
        # It is the jackknife uncertainty of this CV-like REA residual, computed
        # by leaving out one image at a time. This keeps the green line meaningful
        # while adding a real uncertainty estimate for the finite 7-sample ensemble.
        eta_ens = float(sample_std / max(abs(sample_mean), 1.0e-30))
        eta_ens_jack = []
        finite_sample_vals = sample_vals[np.isfinite(sample_vals)]
        n_finite_green = int(finite_sample_vals.size)
        if n_finite_green >= 3:
            for j in range(n_finite_green):
                loo_vals = np.delete(finite_sample_vals, j)
                loo_mean = float(np.nanmean(loo_vals))
                loo_std = float(np.nanstd(loo_vals, ddof=1))
                eta_ens_jack.append(loo_std / max(abs(loo_mean), 1.0e-30))
            eta_ens_jack = np.array(eta_ens_jack, dtype=float)
            eta_ens_jack_mean = float(np.nanmean(eta_ens_jack))
            eta_ens_jack_se = float(
                np.sqrt((n_finite_green - 1) / n_finite_green *
                        np.nansum((eta_ens_jack - eta_ens_jack_mean) ** 2))
            )
        else:
            eta_ens_jack_se = 0.0
        eta_ens_std = eta_ens_jack_se
        eta_ens_sem = eta_ens_jack_se
        eta_C = float(np.nanmean(eta_C_vals)) if not np.all(np.isnan(eta_C_vals)) else np.nan
        eta_C_std = float(np.nanstd(eta_C_vals, ddof=1)) if n_images > 1 and not np.all(np.isnan(eta_C_vals)) else 0.0
        eta_C_sem = float(eta_C_std / np.sqrt(n_images)) if n_images > 1 else 0.0
        # BLUE metric in panel (c): ensemble mean gray-level convergence.
        #
        #     eta_Z(L) = | mean_i[<Z>_L^{(i)}] - mean_i[<Z>_Lref^{(i)}] |
        #                / ( |mean_i[<Z>_Lref^{(i)}]| + eps )
        #
        # This is a REA convergence/bias criterion for the ensemble mean.
        # Its uncertainty is estimated by jackknife leave-one-image-out resampling
        # of the same scalar statistic, avoiding the large and redundant raw
        # image-to-image residual spread.
        delta_Z = abs(sample_mean - ref_ensemble_mean_Z) / max(abs(ref_ensemble_mean_Z), 1.0e-30)
        delta_Z_jack = []
        finite_pair_mask = np.isfinite(sample_vals) & np.isfinite(ref_sample_vals)
        finite_sample_for_blue = sample_vals[finite_pair_mask]
        finite_ref_for_blue = ref_sample_vals[finite_pair_mask]
        n_finite_blue = int(finite_sample_for_blue.size)
        if n_finite_blue >= 3:
            for j in range(n_finite_blue):
                loo_L = np.delete(finite_sample_for_blue, j)
                loo_ref = np.delete(finite_ref_for_blue, j)
                loo_mean_L = float(np.nanmean(loo_L))
                loo_mean_ref = float(np.nanmean(loo_ref))
                delta_Z_jack.append(abs(loo_mean_L - loo_mean_ref) / max(abs(loo_mean_ref), 1.0e-30))
            delta_Z_jack = np.array(delta_Z_jack, dtype=float)
            delta_Z_jack_mean = float(np.nanmean(delta_Z_jack))
            delta_Z_jack_se = float(
                np.sqrt((n_finite_blue - 1) / n_finite_blue *
                        np.nansum((delta_Z_jack - delta_Z_jack_mean) ** 2))
            )
        else:
            delta_Z_jack_se = 0.0
        delta_mean_images = float(np.nanmean(delta_vals)) if not np.all(np.isnan(delta_vals)) else np.nan
        delta_std_images = delta_Z_jack_se
        delta_sem_images = delta_Z_jack_se

        is_reference = L == reference_L

        agg.append({
            "L": L,
            "mean_Z": center_mean,
            "mean_shell": np.nan,
            "sample_mean_Z": sample_mean,
            "sample_std_Z": sample_std,
            "sample_sem_Z": sample_sem,
            "sample_n": n_images,
            # delta_Z_raw is the pointwise ensemble-mean residual. The plotted
            # panel-(c) blue metric delta_Z is replaced after the loop by the
            # 5-point running-median plus monotone tail-maximum plateau of this residual.
            "delta_Z_raw": delta_Z,
            "delta_Z": delta_Z,
            "delta_Z_image_mean": delta_mean_images,
            "delta_Z_std": 0.0,
            "delta_Z_sem": 0.0,
            "eta_C": eta_C,
            "eta_C_std": eta_C_std,
            "eta_C_sem": eta_C_sem,
            "eta_ens": eta_ens,
            "eta_ens_std": eta_ens_std,
            "eta_ens_sem": eta_ens_sem,
            "spectrum_k": rows[0]["spectrum_k"],
            "spectrum_p": rows[0]["spectrum_p"],
            "is_reference": is_reference,
            "image_center_values": center_vals,
            "image_sample_values": sample_vals,
            "image_eta_C_values": eta_C_vals,
            "image_eta_boundary_values": eta_b_vals,
        })

    # BLUE panel-(c) REA metric: running-median plateau envelope.
    # Step 1: compute the raw ensemble-mean residual at each L.
    # Step 2: apply a local 5-point running median to suppress isolated
    # finite-window oscillations without forcing an artificial fit.
    # Step 3: apply a tail-maximum envelope to the smoothed sequence, so
    # the plotted value at L_j is the largest remaining smoothed residual
    # from L_j up to L_ref. This gives a monotone REA convergence criterion:
    # once the curve is below tolerance, all larger tested supports remain
    # below tolerance as well.
    raw = np.array([r.get("delta_Z_raw", r["delta_Z"]) for r in agg], dtype=float)
    is_ref = np.array([bool(r.get("is_reference", False)) for r in agg], dtype=bool)

    # Do not let the reference point participate in the displayed blue
    # convergence curve. The reference residual is zero by construction and
    # otherwise creates an extra blue point where the orange reference residual
    # is hidden.
    raw_for_blue = raw.copy()
    raw_for_blue[is_ref] = np.nan

    med = np.full_like(raw_for_blue, np.nan, dtype=float)
    n_raw = len(raw_for_blue)
    for j in range(n_raw):
        i0 = max(0, j - 2)
        i1 = min(n_raw, j + 3)
        window = raw_for_blue[i0:i1]
        window = window[np.isfinite(window)]
        med[j] = np.nanmedian(window) if window.size else np.nan

    tailmax = np.full_like(med, np.nan, dtype=float)
    running = np.nan
    for j in range(len(med) - 1, -1, -1):
        val = med[j]
        if np.isfinite(val):
            running = val if not np.isfinite(running) else max(running, val)
        tailmax[j] = running

    tailmax[is_ref] = np.nan

    for r, raw_val, med_val, tm in zip(agg, raw, med, tailmax):
        r["delta_Z_raw"] = float(raw_val) if np.isfinite(raw_val) else np.nan
        r["delta_Z_median"] = float(med_val) if np.isfinite(med_val) else np.nan
        r["delta_Z"] = float(tm) if np.isfinite(tm) else np.nan
        r["delta_Z_std"] = 0.0
        r["delta_Z_sem"] = 0.0

    # GREEN panel-(c) REA metric: monotone ensemble-reproducibility envelope.
    # The raw green residual is the coefficient of variation of the seven
    # apparent means at each L.  To use it as a REA decision curve, plot the
    # tail maximum: from L_j onward, what is the largest remaining
    # image-to-image variability?  The green error bar is the jackknife
    # uncertainty of the raw CV at the L value that controls the tail maximum.
    green_raw = np.array([r["eta_ens"] for r in agg], dtype=float)
    green_err_raw = np.array([r["eta_ens_std"] for r in agg], dtype=float)
    green_raw[is_ref] = np.nan
    green_err_raw[is_ref] = np.nan

    green_tail = np.full_like(green_raw, np.nan, dtype=float)
    green_tail_err = np.full_like(green_err_raw, np.nan, dtype=float)
    running_val = np.nan
    running_err = np.nan
    for j in range(len(green_raw) - 1, -1, -1):
        val = green_raw[j]
        err = green_err_raw[j]
        if np.isfinite(val):
            if (not np.isfinite(running_val)) or (val >= running_val):
                running_val = val
                running_err = err if np.isfinite(err) else 0.0
        green_tail[j] = running_val
        green_tail_err[j] = running_err

    green_tail[is_ref] = np.nan
    green_tail_err[is_ref] = np.nan
    for r, gv, ge in zip(agg, green_tail, green_tail_err):
        r["eta_ens_raw"] = r["eta_ens"]
        r["eta_ens"] = float(gv) if np.isfinite(gv) else np.nan
        r["eta_ens_std"] = float(ge) if np.isfinite(ge) else np.nan
        r["eta_ens_sem"] = float(ge) if np.isfinite(ge) else np.nan

    return agg


def build_ensemble_spectra(per_image_results: list[list[dict]], valid_window_sizes: list[int]):
    """Build ensemble-averaged spectra for panel (d).

    For each L, spectra from all input images are interpolated to the k-grid
    of the first image and averaged. The returned error field is the
    image-to-image standard deviation or standard error, depending on
    --errorbar-mode.
    """
    n_images = len(per_image_results)
    spectra = []
    for idx, L in enumerate(valid_window_sizes):
        k0 = np.asarray(per_image_results[0][idx]["spectrum_k"], dtype=float)
        p_stack = []
        for image_results in per_image_results:
            k = np.asarray(image_results[idx]["spectrum_k"], dtype=float)
            p = np.asarray(image_results[idx]["spectrum_p"], dtype=float)
            mask = np.isfinite(k) & np.isfinite(p) & (k > 0) & (p > 0)
            if np.count_nonzero(mask) < 3:
                p_interp = np.full_like(k0, np.nan, dtype=float)
            else:
                # Interpolate in log-amplitude space to keep positive spectra stable.
                p_interp = np.exp(np.interp(k0, k[mask], np.log(p[mask]), left=np.nan, right=np.nan))
            p_stack.append(p_interp)
        p_stack = np.vstack(p_stack)
        p_mean = np.nanmean(p_stack, axis=0)
        p_std = np.nanstd(p_stack, axis=0, ddof=1) if n_images > 1 else np.zeros_like(p_mean)
        p_sem = p_std / np.sqrt(n_images) if n_images > 1 else np.zeros_like(p_mean)
        spectra.append({
            "L": L,
            "spectrum_k": k0,
            "spectrum_p_mean": p_mean,
            "spectrum_p_std": p_std,
            "spectrum_p_sem": p_sem,
        })
    return spectra


def make_combined_figure(z_representative, ensemble_results, spectra_results, accepted_L, no_pdf=False):
    """Make the final manuscript-style 4-panel figure plus a separate
    consistency-check figure.

    Main figure:
      (a) representative BSE/gray field with selected REA square and colorbar;
      (b) ensemble mean \bar{Z}_L only;
      (c) primary structural REA metric only, eta_C/tau_C;
      (d) second-order spectral collapse.

    Separate consistency figure:
      (a) mean gray-level residual eta_Z/tau_Z;
      (b) ensemble reproducibility residual eta_ens/tau_ens.
    """
    Lvals = np.array([r["L"] for r in ensemble_results], dtype=float)
    sample_mean_Z = np.array([r["sample_mean_Z"] for r in ensemble_results], dtype=float)

    raw_delta = diagnostic_array(ensemble_results, "delta_Z", hide_reference=True)
    raw_etaC = diagnostic_array(ensemble_results, "eta_C", hide_reference=False)
    raw_etaens = diagnostic_array(ensemble_results, "eta_ens", hide_reference=False)

    raw_etaC_err = diagnostic_array(
        ensemble_results,
        "eta_C_sem" if ERRORBAR_MODE == "sem" else "eta_C_std",
        hide_reference=True,
    )
    raw_etaens_err = diagnostic_array(
        ensemble_results,
        "eta_ens_sem" if ERRORBAR_MODE == "sem" else "eta_ens_std",
        hide_reference=False,
    )

    crit_delta = raw_delta / TAU_Z
    crit_etaC = raw_etaC / TAU_C
    crit_etaens = raw_etaens / TAU_ENS
    crit_etaC_err = raw_etaC_err / TAU_C
    crit_etaens_err = raw_etaens_err / TAU_ENS

    # The reference spectral residual is zero by construction. To make the
    # reference support visible on a log-scale panel, plot zero positive values
    # at a very small display floor without changing the numerical output/CSV.
    crit_etaC_plot = np.array(crit_etaC, dtype=float, copy=True)
    zero_etaC = np.isfinite(crit_etaC_plot) & (crit_etaC_plot <= 0.0)
    crit_etaC_plot[zero_etaC] = LOG_FLOOR

    fig, axes = plt.subplots(2, 2, figsize=(7.85, 6.35))

    # --------------------------------------------------------
    # Panel (a): representative image and selected REA window
    # --------------------------------------------------------
    ax = axes[0, 0]
    im = ax.imshow(z_representative, cmap="gray", origin="lower")
    box_L = accepted_L if accepted_L is not None else int(Lvals[-1])
    x0, y0, w, h = centered_window_bounds(z_representative.shape, int(box_L))
    rect = plt.Rectangle((x0, y0), w, h, fill=False, edgecolor="white", linewidth=1.5)
    ax.add_patch(rect)
    ax.text(
        0.03, 0.95, "(a)", transform=ax.transAxes, ha="left", va="top",
        fontsize=11, fontweight="bold", color="white",
        bbox=dict(facecolor="black", alpha=0.35, edgecolor="none", pad=2.0),
    )
    ax.set_xlabel("x [pixels]")
    ax.set_ylabel("y [pixels]")
    # Keep the gray-level colorbar, but place it far enough from the axes so it
    # does not overlap the x-axis tick labels in the compact 4-panel layout.
    cbar = fig.colorbar(
        im,
        ax=ax,
        orientation="vertical",
        fraction=0.035,
        pad=0.070,
        shrink=0.74,
    )
    cbar.set_label(r"$Z_{\rm gray}$", fontsize=8)
    cbar.ax.tick_params(labelsize=7, pad=1.5)

    # --------------------------------------------------------
    # Panel (b): ensemble mean of accepted stationary windows.
    #
    # The raw curve can be noisy because each L samples a finite collection of
    # stationary windows. For the manuscript panel, an optional running-median
    # trend and a reference tolerance band make the plateau region easier to
    # read without changing the REA decision, which remains panel (c).
    # --------------------------------------------------------
    ax = axes[0, 1]
    panel_b_smooth = running_nanmedian_1d(sample_mean_Z, PANEL_B_SMOOTH_WINDOW)

    if PANEL_B_REFERENCE_BAND:
        # Visual band: finite reference mean +/- tau_Z.  This is a display guide
        # for the apparent mean, not the primary REA criterion.
        if accepted_L is not None and np.any(Lvals == accepted_L):
            ref_for_band = float(sample_mean_Z[np.where(Lvals == accepted_L)[0][0]])
        else:
            finite_mean = sample_mean_Z[np.isfinite(sample_mean_Z)]
            ref_for_band = float(finite_mean[-1]) if finite_mean.size else np.nan
        if np.isfinite(ref_for_band):
            ax.axhspan(
                ref_for_band * (1.0 - TAU_Z),
                ref_for_band * (1.0 + TAU_Z),
                color="C0",
                alpha=0.08,
                linewidth=0.0,
                label=rf"$\pm\tau_Z$ band",
            )

    if PANEL_B_MODE in {"raw", "raw+smooth"}:
        ax.plot(
            Lvals,
            sample_mean_Z,
            "o-",
            linewidth=1.15 if PANEL_B_MODE == "raw+smooth" else 1.45,
            markersize=3.5,
            label=r"$\overline{Z}_L$",
            color="C0",
            alpha=0.35 if PANEL_B_MODE == "raw+smooth" else 0.98,
        )

    if PANEL_B_MODE in {"smooth", "raw+smooth"}:
        ax.plot(
            Lvals,
            panel_b_smooth,
            "o-",
            linewidth=1.75,
            markersize=3.8,
            label="median trend",
            color="C0",
            alpha=1.0,
        )

    add_rea_line(ax, accepted_L)
    ax.text(
        0.03, 0.95, "(b)", transform=ax.transAxes, ha="left", va="top",
        fontsize=11, fontweight="bold",
    )
    ax.set_xlabel(r"window size $L$ [pixels]")
    ax.set_ylabel(r"ensemble mean $\overline{Z}_L$")
    ax.grid(True, alpha=0.22)
    ax.legend(frameon=False, fontsize=7.0, loc="best")

    # --------------------------------------------------------
    # Panel (c): primary structural metric only
    # --------------------------------------------------------
    ax = axes[1, 0]
    if accepted_L is not None:
        err_mask = Lvals >= accepted_L
    else:
        err_mask = np.ones_like(Lvals, dtype=bool)

    finite_eta = crit_etaC_plot[np.isfinite(crit_etaC_plot) & (crit_etaC_plot > 0)]
    if finite_eta.size:
        panel_c_ymin = max(min(np.nanmin(finite_eta) / 2.0, 0.8), LOG_FLOOR)
        panel_c_ymax = max(np.nanmax(finite_eta) * 1.8, 1.35)
    else:
        panel_c_ymin = LOG_FLOOR
        panel_c_ymax = 1.5

    ax.fill_between(
        [np.nanmin(Lvals), np.nanmax(Lvals)],
        panel_c_ymin,
        1.0,
        color="0.92",
        alpha=0.70,
        zorder=0,
    )
    # Optional display smoothing for the structural residual. The REA decision
    # and CSV output remain based on the unsmoothed eta_C values.  The raw
    # values can be shown faintly, while the running-median log-space trend
    # makes the post-REA plateau easier to read.
    panel_c_smooth = running_log_median_positive(crit_etaC_plot, PANEL_C_SMOOTH_WINDOW)

    if PANEL_C_MODE in {"raw", "raw+smooth"}:
        ax.plot(
            Lvals,
            crit_etaC_plot,
            "s-",
            linewidth=1.05 if PANEL_C_MODE == "raw+smooth" else 1.55,
            markersize=3.6,
            label=r"$\eta_{\widehat C}/\tau_{\widehat C}$",
            color="C1",
            alpha=0.35 if PANEL_C_MODE == "raw+smooth" else 1.0,
        )

    if PANEL_C_MODE in {"smooth", "raw+smooth"}:
        ax.plot(
            Lvals,
            panel_c_smooth,
            "s-",
            linewidth=1.80,
            markersize=4.0,
            label=rf"$\eta_{{\widehat C}}/\tau_{{\widehat C}}$ median trend",
            color="C1",
            alpha=1.0,
        )

    # Error bars are shown only from L_REA onward, where the reader needs to
    # assess the post-selection stability. This keeps the pre-REA transient
    # from visually dominating panel (c).
    err_mask_orange = err_mask & np.isfinite(crit_etaC_plot) & np.isfinite(crit_etaC_err)
    ax.errorbar(
        Lvals[err_mask_orange],
        crit_etaC_plot[err_mask_orange],
        yerr=crit_etaC_err[err_mask_orange],
        fmt="none",
        ecolor="C1",
        elinewidth=1.0,
        capsize=2.5,
        capthick=1.0,
        zorder=2.5,
        alpha=0.90,
    )
    ax.set_yscale("log")
    ax.axhline(1.0, color="black", linestyle=":", linewidth=1.1)
    add_rea_line(ax, accepted_L)
    ax.text(
        0.03, 0.95, "(c)", transform=ax.transAxes, ha="left", va="top",
        fontsize=11, fontweight="bold",
    )
    ax.set_ylim(panel_c_ymin, panel_c_ymax)
    ax.set_xlabel(r"window size $L$ [pixels]")
    ax.set_ylabel("structural residual / tolerance")
    if np.any(np.isfinite(Lvals)):
        xmin = float(np.nanmin(Lvals))
        xmax = float(np.nanmax(Lvals))
        pad = 0.07 * max(xmax - xmin, 1.0)
        ax.set_xlim(xmin - 0.02 * max(xmax - xmin, 1.0), xmax + pad)
    ax.grid(True, which="major", alpha=0.24)
    ax.grid(True, which="minor", alpha=0.10)
    ax.legend(frameon=False, fontsize=7.3, handlelength=2.0)

    # --------------------------------------------------------
    # Panel (d): spectral collapse
    # --------------------------------------------------------
    ax = axes[1, 1]
    d_xvals = []
    d_yvals = []
    shown = set(SELECTED_L_FOR_SPECTRA)
    for r in spectra_results:
        if r["L"] in shown:
            k = np.asarray(r["spectrum_k"], dtype=float)
            p = np.asarray(r.get("spectrum_p_mean", r.get("spectrum_p")), dtype=float)
            perr = np.asarray(
                r.get("spectrum_p_sem" if ERRORBAR_MODE == "sem" else "spectrum_p_std", np.zeros_like(p)),
                dtype=float,
            )
            mask = (k > 0) & (p > 0) & np.isfinite(k) & np.isfinite(p)
            if np.any(mask):
                norm = np.nanmax(p[mask]) + 1.0e-30
                y = p[mask] / norm
                yerr_spec = perr[mask] / norm
                d_xvals.append(k[mask])
                d_yvals.append(y)
                ax.loglog(k[mask], y, linewidth=1.25, label=rf"$L={r['L']}$")
                lo = np.maximum(y - yerr_spec, 1.0e-30)
                hi = y + yerr_spec
                d_yvals.append(lo)
                d_yvals.append(hi)
                ax.fill_between(k[mask], lo, hi, alpha=0.10, linewidth=0)

    if d_xvals and d_yvals:
        dx = np.concatenate([np.asarray(v, dtype=float) for v in d_xvals])
        dy = np.concatenate([np.asarray(v, dtype=float) for v in d_yvals])
        dx = dx[np.isfinite(dx) & (dx > 0)]
        dy = dy[np.isfinite(dy) & (dy > 0)]
        if dx.size:
            ax.set_xlim(dx.min() / 1.15, dx.max() * 1.15)
        if dy.size:
            ax.set_ylim(max(dy.min() / 1.8, 1.0e-8), min(dy.max() * 1.45, 2.0))
    ax.text(
        0.03, 0.95, "(d)", transform=ax.transAxes, ha="left", va="top",
        fontsize=11, fontweight="bold",
    )
    ax.set_xlabel(r"$\log_{10}(k\,[\mathrm{pixel}^{-1}])$")
    ax.set_ylabel(r"$\log_{10}\!\left(\widehat{C}_L(k)\right)$")
    ax.grid(True, which="both", alpha=0.22)
    ax.legend(frameon=False, fontsize=6.6, ncol=2, handlelength=1.8, columnspacing=0.9)

    fig.tight_layout(pad=0.9, w_pad=1.0, h_pad=0.9)
    save_png_pdf(
        fig,
        f"{OUTPUT_PREFIX}_combined_figure4_style_summary",
        dpi=SAVE_DPI,
        save_pdf=not no_pdf,
    )
    plt.close(fig)

    # --------------------------------------------------------
    # Separate consistency-check figure: blue and green diagnostics.
    # This panel is intentionally kept in the v20 scientific style:
    #   - no smoothing;
    #   - no error bars;
    #   - no panel-b/panel-c display transformations;
    #   - only the two auxiliary diagnostics used in the manuscript.
    # These are consistency checks, not primary REA criteria.
    # --------------------------------------------------------
    fig2, ax = plt.subplots(1, 1, figsize=(4.45, 3.15))

    consistency_vals = np.concatenate([
        crit_delta[np.isfinite(crit_delta) & (crit_delta > 0)],
        crit_etaens[np.isfinite(crit_etaens) & (crit_etaens > 0)],
    ])
    if consistency_vals.size:
        # Use a stable lower display cap so the v20-style plot does not collapse
        # vertically when all auxiliary residuals are well below one.
        ymin2 = max(min(np.nanmin(consistency_vals) / 2.0, 8.0e-2), LOG_FLOOR)
        ymax2 = max(np.nanmax(consistency_vals) * 1.8, 1.35)
    else:
        ymin2, ymax2 = LOG_FLOOR, 1.5

    ax.fill_between(
        [np.nanmin(Lvals), np.nanmax(Lvals)],
        ymin2,
        1.0,
        color="0.92",
        alpha=0.70,
        zorder=0,
    )
    ax.plot(
        Lvals,
        crit_delta,
        "o-",
        linewidth=1.45,
        markersize=4.0,
        color="C0",
        label=r"$\eta_Z^{\rm med\text{-}tail}/\tau_Z$",
    )
    ax.plot(
        Lvals,
        crit_etaens,
        "^-",
        linewidth=1.45,
        markersize=4.0,
        color="C2",
        label=r"$\eta_{\rm ens}^{\rm tail}/\tau_{\rm ens}$",
    )
    ax.set_yscale("log")
    ax.axhline(1.0, color="black", linestyle=":", linewidth=1.0)
    add_rea_line(ax, accepted_L)
    ax.set_xlabel(r"window size $L$ [pixels]")
    ax.set_ylabel("consistency residual / tolerance")
    ax.set_ylim(ymin2, ymax2)
    ax.grid(True, which="major", alpha=0.24)
    ax.grid(True, which="minor", alpha=0.10)
    ax.legend(frameon=False, fontsize=7.0, loc="best")

    fig2.tight_layout(pad=0.9)
    save_png_pdf(
        fig2,
        f"{OUTPUT_PREFIX}_consistency_checks",
        dpi=SAVE_DPI,
        save_pdf=not no_pdf,
    )
    plt.close(fig2)

def save_combined_csv(path: str, image_names: list[str], ensemble_results: list[dict]):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("L,ensemble_center_mean_Z,ensemble_sample_mean_Z,ensemble_sample_std_Z,ensemble_sample_sem_Z,n_images,eta_Z,eta_C,eta_ens,is_reference")
        for name in image_names:
            stem = Path(name).stem.replace(",", "_")
            fh.write(f",sample_mean_Z__{stem}")
        fh.write("\n")
        for r in ensemble_results:
            fh.write(
                f"{r['L']},{r['mean_Z']:.12e},{r['sample_mean_Z']:.12e},"
                f"{r['sample_std_Z']:.12e},{r['sample_sem_Z']:.12e},{r['sample_n']},"
                f"{r['delta_Z']:.12e},{r['eta_C']:.12e},{r['eta_ens']:.12e},{r['is_reference']}"
            )
            for val in r["image_sample_values"]:
                fh.write(f",{val:.12e}")
            fh.write("\n")


def run_combined_images(args):
    """Process all input images as one ensemble and save one combined figure."""
    global OUTPUT_PREFIX
    prefix = args.output_prefix if args.output_prefix is not None else "combined_bse_REA"
    OUTPUT_PREFIX = str(OUTPUT_DIR / prefix)

    input_paths = [Path(p) for p in args.images]
    for p in input_paths:
        if not p.exists():
            raise FileNotFoundError(f"Could not find input image: {p.resolve()}")

    print("\nCombined ensemble REA mode")
    print("=" * 72)
    print(f"Number of images             : {len(input_paths)}")
    print("Images:")
    for i, p in enumerate(input_paths, start=1):
        print(f"  {i}: {p}")
    print("=" * 72)

    gray_fields = []
    z_fields = []
    masks = []
    dx_values = []
    shapes = []
    for p in input_paths:
        gray, dx_eff, _, _ = load_gray_image_as_array(p, crop_box=CROP_BOX, max_side=MAX_SIDE_FOR_ANALYSIS)
        z = gray_to_analysis_field(gray)
        stationary_mask = load_stationary_mask_for_image(p, target_shape=z.shape)
        gray_fields.append(gray)
        z_fields.append(z)
        masks.append(stationary_mask)
        dx_values.append(dx_eff)
        shapes.append(z.shape)

    if np.nanmax(dx_values) - np.nanmin(dx_values) > 1e-12:
        print("Warning: effective pixel sizes are not identical across images. Using the first image for spectral units.")
    dx_eff = dx_values[0]

    valid_window_sizes = common_valid_window_sizes_with_masks(shapes, masks)
    reference_L = REFERENCE_L if REFERENCE_L in valid_window_sizes else valid_window_sizes[-1]
    if reference_L != REFERENCE_L:
        print(f"Warning: requested reference-L={REFERENCE_L} is not common to all images. Using {reference_L}.")

    print("\nCommon analysis setup")
    print("=" * 72)
    print(f"Common window sizes          : {valid_window_sizes}")
    print(f"Reference window             : {reference_L} pixels")
    print(f"Field mode                   : {FIELD_MODE}")
    print(f"Stationary masks             : {'enabled' if STATIONARY_MASK_DIR is not None else 'disabled'}")
    if STATIONARY_MASK_DIR is not None:
        print(f"Mask directory               : {STATIONARY_MASK_DIR}")
        print(f"Minimum mask fraction/window : {MIN_MASK_FRACTION}")
    print(f"Error bars                   : image-to-image {ERRORBAR_MODE}")
    print(f"Panel-(c) Z metric           : 5-point running-median + monotone tail-max plateau, reference point hidden of ensemble mean gray-level residual, line only")
    print(f"Panel-(c) spectral metric    : {SPECTRAL_BAND}-k relative residual, weight={LOW_K_WEIGHT}")
    print("=" * 72)

    per_image_results = []
    for i, (z, m) in enumerate(zip(z_fields, masks), start=1):
        print(f"\nComputing image {i}/{len(z_fields)}")
        results, _ = compute_image_results(z, dx_eff, valid_window_sizes, reference_L, stationary_mask=m)
        per_image_results.append(results)

    ensemble_results = aggregate_ensemble_results(per_image_results, valid_window_sizes, reference_L)
    accepted_L = find_persistent_rea(
        ensemble_results,
        tau_Z=TAU_Z,
        tau_C=TAU_C,
        tail_fraction=TAIL_FRACTION,
    )

    print("\nCombined ensemble table")
    print("=" * 120)
    print(f"{'L':>6s} {'<Z>_ens':>12s} {'err':>12s} {'eta_Z':>12s} {'eta_C':>12s} {'eta_ens':>12s} {'ok':>8s} {'reference':>10s}")
    print("-" * 120)
    for r in ensemble_results:
        err = r["sample_sem_Z"] if ERRORBAR_MODE == "sem" else r["sample_std_Z"]
        ok = check_window_ok(r, TAU_Z, TAU_C)
        print(
            f"{r['L']:6d} {r['sample_mean_Z']:12.5e} {err:12.4e} "
            f"{r['delta_Z']:12.4e} {r['eta_C'] if not np.isnan(r['eta_C']) else np.nan:12.4e} "
            f"{r['eta_ens']:12.4e} {str(ok):>8s} {str(r['is_reference']):>10s}"
        )
    print("=" * 120)

    if accepted_L is None:
        print("No robust persistent ensemble L_REA found with the selected tolerances.")
    else:
        print(f"Robust persistent ensemble L_REA = {accepted_L} pixels")

    rep_idx = max(0, min(int(args.representative_image) - 1, len(z_fields) - 1))
    spectra_results = build_ensemble_spectra(per_image_results, valid_window_sizes)
    make_combined_figure(z_fields[rep_idx], ensemble_results, spectra_results, accepted_L, no_pdf=args.no_pdf)

    csv_path = f"{OUTPUT_PREFIX}_combined_window_statistics.csv"
    save_combined_csv(csv_path, [str(p) for p in input_paths], ensemble_results)
    print("\nSaved combined outputs:")
    if args.no_pdf:
        print(f"  {OUTPUT_PREFIX}_combined_figure4_style_summary.png")
    else:
        print(f"  {OUTPUT_PREFIX}_combined_figure4_style_summary.png / .pdf")
    print(f"  {csv_path}")


# ============================================================
# Main
# ============================================================


def output_prefix_for_image(image_path: Path, args, n_images: int) -> str:
    safe_stem = image_path.stem.replace(" ", "_")
    if args.output_prefix is None:
        prefix = safe_stem
    elif n_images == 1:
        prefix = args.output_prefix
    else:
        prefix = f"{args.output_prefix}_{safe_stem}"
    return str(OUTPUT_DIR / prefix)


def run_one_image(input_path: Path, output_prefix: str, args=None):
    global OUTPUT_PREFIX
    OUTPUT_PREFIX = output_prefix

    if not input_path.exists():
        raise FileNotFoundError(f"Could not find input image: {input_path.resolve()}")

    gray, dx_eff, _, _ = load_gray_image_as_array(
        input_path,
        crop_box=CROP_BOX,
        max_side=MAX_SIDE_FOR_ANALYSIS,
    )
    z_full = gray_to_analysis_field(gray)

    valid_window_sizes = validate_window_sizes(z_full.shape)
    reference_L = REFERENCE_L if REFERENCE_L in valid_window_sizes else valid_window_sizes[-1]

    kmax = np.pi / dx_eff
    q_grid = np.linspace(0.0, kmax, N_COMMON_K_GRID)
    weights = trapezoid_weights(q_grid)

    # Fixed sample-level scale used only for the green boundary-shell metric.
    # This prevents the green residual from being dominated by an L-dependent
    # denominator while preserving the boundary-vs-full-window numerator.
    boundary_material_scale = robust_material_scale(z_full)

    spectra = {}
    for L in valid_window_sizes:
        z_L = centered_square_window(z_full, L)
        spectra[L] = radial_average_spectrum(z_L, dx=dx_eff)
    k_ref, p_ref = spectra[reference_L]
    ref_sample_mean_Z, _, _, _ = window_mean_statistics(z_full, reference_L, approx_n=N_WINDOW_SAMPLES)
    ref_center_field = centered_square_window(z_full, reference_L)
    _, ref_center_mean_Z, _ = boundary_shell_diagnostic(
        ref_center_field, alpha=ALPHA_SHELL, material_scale=boundary_material_scale
    )

    results = []
    previous_mean_Z = None

    print("\nBSE gray-level REA scan")
    print("=" * 118)
    print(f"Output prefix          = {OUTPUT_PREFIX}")
    print(f"Field mode             = {FIELD_MODE}")
    print(f"k_max = pi / Delta x_eff = {kmax:.6e}")
    print(f"Reference window L_ref = {reference_L} pixels")
    print("=" * 118)
    print(
        f"{'L':>6s} "
        f"{'mean_Z':>12s} {'avg_Z':>12s} {'err_Z':>12s} "
        f"{'eta_Z':>12s} {'eta_C':>12s} {'eta_shell':>12s} "
        f"{'ok':>8s} {'reference':>10s}"
    )
    print("-" * 118)

    for L in valid_window_sizes:
        z_L = centered_square_window(z_full, L)
        k, p = spectra[L]
        eta_bnd, mean_Z, mean_shell = boundary_shell_diagnostic(
            z_L, alpha=ALPHA_SHELL, material_scale=boundary_material_scale
        )
        sample_mean_Z, sample_std_Z, sample_sem_Z, sample_n = window_mean_statistics(
            z_full, L, approx_n=N_WINDOW_SAMPLES
        )
        sample_err_Z = sample_sem_Z if ERRORBAR_MODE == "sem" else sample_std_Z

        is_reference = L == reference_L
        if is_reference:
            eta_C = np.nan
        else:
            eta_C = low_k_spectral_residual_full_band(k, p, k_ref, p_ref, q_grid, weights)

        delta_Z = abs(sample_mean_Z - ref_sample_mean_Z) / max(abs(ref_sample_mean_Z), 1.0e-30)
        eta_Z_centered = abs(mean_Z - ref_center_mean_Z) / max(abs(ref_center_mean_Z), 1.0e-30)

        r = {
            "L": L,
            "mean_Z": mean_Z,
            "mean_shell": mean_shell,
            "sample_mean_Z": sample_mean_Z,
            "sample_std_Z": sample_std_Z,
            "sample_sem_Z": sample_sem_Z,
            "sample_n": sample_n,
            "delta_Z": delta_Z,
            "eta_Z": delta_Z,
            "eta_Z_centered": eta_Z_centered,
            "eta_C": eta_C,
            "eta_boundary": eta_bnd,
            "spectrum_k": k,
            "spectrum_p": p,
            "is_reference": is_reference,
        }
        results.append(r)
        ok = check_window_ok(r, TAU_Z, TAU_C)

        print(
            f"{L:6d} "
            f"{mean_Z:12.5e} {sample_mean_Z:12.5e} {sample_err_Z:12.4e} "
            f"{delta_Z:12.4e} "
            f"{eta_C if not np.isnan(eta_C) else np.nan:12.4e} "
            f"{eta_bnd:12.4e} "
            f"{str(ok):>8s} {str(is_reference):>10s}"
        )
        previous_mean_Z = mean_Z

    print("-" * 118)

    print("\nNormalized convergence diagnostics")
    print("=" * 95)
    print(
        f"{'L':>6s} "
        f"{'deltaZ/tau':>14s} "
        f"{'etaC/tau':>14s} "
        f"{'etab/tau':>14s} "
        f"{'all below 1':>14s}"
    )
    print("-" * 95)
    for r in results:
        if r["is_reference"] or np.isnan(r["delta_Z"]) or np.isnan(r["eta_C"]):
            continue
        d1 = r["delta_Z"] / TAU_Z
        d2 = r["eta_C"] / TAU_C
        d3 = r["eta_boundary"] / TAU_SHELL
        ok_all = d1 <= 1.0 and d2 <= 1.0 and d3 <= 1.0
        print(f"{r['L']:6d} {d1:14.4e} {d2:14.4e} {d3:14.4e} {str(ok_all):>14s}")
    print("=" * 95)

    accepted_L = find_persistent_rea(
        results,
        tau_Z=TAU_Z,
        tau_C=TAU_C,
        tail_fraction=TAIL_FRACTION,
    )

    if accepted_L is None:
        print("No robust persistent L_REA found with the selected tolerances.")
        print("The tested windows do not prove REA convergence for this image/tolerance set.")
    else:
        print(f"Robust persistent L_REA = {accepted_L} pixels")
        print(f"REA side length         = {accepted_L * dx_eff:.6g} in the supplied pixel-size units")


    # Save numerical values used in panels, including average values and error bars.
    csv_path = f"{OUTPUT_PREFIX}_window_statistics.csv"
    with open(csv_path, "w", encoding="utf-8") as fh:
        fh.write("L,center_mean_Z,sample_mean_Z,sample_std_Z,sample_sem_Z,sample_n,delta_Z,eta_C,eta_boundary,is_reference\n")
        for r in results:
            fh.write(
                f"{r['L']},{r['mean_Z']:.12e},{r['sample_mean_Z']:.12e},"
                f"{r['sample_std_Z']:.12e},{r['sample_sem_Z']:.12e},{r['sample_n']},"
                f"{r['delta_Z']:.12e},{r['eta_C']:.12e},{r['eta_boundary']:.12e},{r['is_reference']}\n"
            )
    print(f"Saved statistics table: {csv_path}")

    make_figures(z_full=z_full, gray_full=gray, results=results, accepted_L=accepted_L, kmax=kmax, save_diagnostics=args.save_diagnostics, no_pdf=args.no_pdf)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    apply_args_to_globals(args)

    # Manuscript default: all input images are treated as independent fields of view and produce one figure.
    # Per-image figures are produced only if explicitly requested with --separate-images.
    if not args.separate_images:
        run_combined_images(args)
        return

    n_images = len(args.images)
    for image_name in args.images:
        input_path = Path(image_name)
        out_prefix = output_prefix_for_image(input_path, args, n_images)
        run_one_image(input_path, out_prefix, args=args)


if __name__ == "__main__":
    main()
