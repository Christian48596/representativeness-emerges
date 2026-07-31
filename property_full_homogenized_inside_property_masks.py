#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
property_full_homogenized_inside_property_masks.py
=================================================

Full property-level homogenization inside stationary masks defined on the
property-map grid.

This script is intended to reproduce the *type* of apparent-property curves
computed in the MATLAB workflow:

    kappa_app_fem.m              -> kappa_x^app(L), kappa_y^app(L)
    stiffness_tensor.m+c_app_fem.m -> C_app(L), E_x^app(L), E_y^app(L)

but with candidate windows restricted by a stationary mask that was computed
on the same property-map grid, e.g. with detect_stationary_domains_general.py.

Inputs
------
CSV/NPY property maps:
    kappa(x)  [W m^-1 K^-1]
    E(x)      [GPa]
    nu(x)     [-]

Mask images/arrays:
    S(x) in {0,1}, same grid as the property maps.

Outputs
-------
    *_kappa_app_curves.png/pdf
    *_C_app_curves.png/pdf
    *_E_app_curves.png/pdf
    *_summary.csv
    *_window_values.csv

Scientific notes
----------------
1. Use numerical CSV/NPY property maps, not TIFF/PNG visualization maps. The
   MATLAB pipeline may export TIFFs after mat2gray()/uint16 scaling; those do
   not contain physical values.
2. This script implements a periodic Q4 scalar-conduction solver for kappa and
   a periodic Q4 vector-elasticity solver for C_app. The elasticity part builds
   the local isotropic stiffness from E and nu using plane strain or plane stress.
3. The apparent Young moduli are derived from the compliance matrix:
       S_app = inv(C_app),  E_x=1/S_11, E_y=1/S_22.
4. Full elasticity is much more expensive than scalar averaging. Start with a
   small number of samples and modest maximum L, then increase.
5. The value supplied through --L-rea is the structural representative
   elementary area transferred to the QEMSCAN property-map grid. It is plotted
   as an independently predicted support and is not fitted to the apparent
   thermal or elastic property curves.

Example
-------
python property_full_homogenized_inside_property_masks.py \
  --kappa-maps kappaMap_1.csv ... kappaMap_7.csv \
  --E-maps EMap_1.csv ... EMap_7.csv \
  --nu-maps nuMap_1.csv ... nuMap_7.csv \
  --mask-dir stationary_property_combined_STRICT \
  --mask-template 'property_combined_{i}_stationary_mask.tif' \
  --window-sizes 32,48,64,80,96,112,128,160,192,204,224,256,320,384,448,512 \
  --L-rea 204 \
  --property-pixel-size-mm 0.01 \
  --n-window-samples 4 \
  --min-mask-fraction 0.98 \
  --elastic-mode plane_strain
"""
from __future__ import annotations

import argparse
import csv
import math
import time
from pathlib import Path
from typing import List, Sequence, Tuple, Dict, Optional

import numpy as np
from PIL import Image
Image.MAX_IMAGE_PIXELS = None

import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as spla

EPS = 1e-30

# =============================================================================
# I/O
# =============================================================================

def load_numeric_map(path: str | Path) -> np.ndarray:
    """Load a 2D numeric array from csv/txt/dat/npy/npz or image."""
    path = Path(path)
    suf = path.suffix.lower()
    if suf == ".npy":
        arr = np.load(path)
    elif suf == ".npz":
        d = np.load(path)
        if len(d.files) != 1:
            raise ValueError(f"{path} contains multiple arrays: {d.files}")
        arr = d[d.files[0]]
    elif suf in {".csv", ".txt", ".dat"}:
        arr = np.genfromtxt(path, delimiter=",")
        if arr.ndim != 2 or np.isnan(arr).all():
            arr = np.genfromtxt(path)
    else:
        arr = np.asarray(Image.open(path))
        if arr.ndim == 3:
            arr = arr[..., :3].astype(float).mean(axis=2)
    arr = np.asarray(arr, dtype=np.float64)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D map for {path}; got {arr.shape}")
    if not np.isfinite(arr).all():
        raise ValueError(f"Map {path} contains non-finite values")
    return arr


def load_mask(path: str | Path) -> np.ndarray:
    """Load binary mask and return bool array."""
    return load_numeric_map(path) > 0


def parse_int_list(text: str) -> List[int]:
    vals = [int(x.strip()) for x in text.split(',') if x.strip()]
    if not vals:
        raise ValueError("empty --window-sizes list")
    return sorted(set(vals))


def resolve_masks(args, n: int, stems: Sequence[str]) -> List[Path]:
    """Resolve mask filenames from explicit --masks or --mask-dir/template."""
    if args.masks is not None:
        if len(args.masks) != n:
            raise ValueError("--masks must have same length as --kappa-maps")
        return [Path(p) for p in args.masks]
    if args.mask_dir is None:
        raise ValueError("Provide --masks or --mask-dir plus --mask-template")
    out = []
    for i, stem in enumerate(stems, start=1):
        out.append(Path(args.mask_dir) / args.mask_template.format(i=i, stem=stem))
    return out

# =============================================================================
# Sampling inside stationary masks
# =============================================================================

def integral_image(arr: np.ndarray) -> np.ndarray:
    return np.pad(arr.astype(float).cumsum(0).cumsum(1), ((1, 0), (1, 0)), mode="constant")


def rect_sum(ii: np.ndarray, y0: int, x0: int, L: int) -> float:
    return float(ii[y0+L, x0+L] - ii[y0, x0+L] - ii[y0+L, x0] + ii[y0, x0])


def sample_valid_windows(mask: np.ndarray, L: int, n_samples: int, min_fraction: float,
                         rng: np.random.Generator, max_random_tries: int = 100000) -> List[Tuple[int, int, float]]:
    """Sample square LxL windows with mask coverage >= min_fraction."""
    H, W = mask.shape
    if L > H or L > W:
        return []
    ii = integral_image(mask)
    area = float(L * L)
    max_y, max_x = H - L, W - L
    accepted: List[Tuple[int, int, float]] = []
    seen = set()
    tries = 0
    while len(accepted) < n_samples and tries < max_random_tries:
        tries += 1
        y0 = int(rng.integers(0, max_y + 1))
        x0 = int(rng.integers(0, max_x + 1))
        if (y0, x0) in seen:
            continue
        seen.add((y0, x0))
        frac = rect_sum(ii, y0, x0, L) / area
        if frac >= min_fraction:
            accepted.append((y0, x0, frac))
    # deterministic fallback on a grid
    if len(accepted) < n_samples:
        step = max(1, L // 4)
        for y0 in range(0, max_y + 1, step):
            for x0 in range(0, max_x + 1, step):
                if (y0, x0) in seen:
                    continue
                frac = rect_sum(ii, y0, x0, L) / area
                if frac >= min_fraction:
                    accepted.append((y0, x0, frac))
                    if len(accepted) >= n_samples:
                        return accepted
    return accepted

# =============================================================================
# Periodic Q4 scalar conductivity solver
# =============================================================================

def q4_scalar_local(hx: float, hy: float):
    g = np.array([-1.0, 1.0]) / math.sqrt(3.0)
    Ke0 = np.zeros((4, 4), float)
    fe0x = np.zeros(4, float)
    fe0y = np.zeros(4, float)
    dNdx_list = []
    dNdy_list = []
    weights = []
    for xi in g:
        for eta in g:
            dN_dxi = 0.25 * np.array([-(1-eta), (1-eta), (1+eta), -(1+eta)])
            dN_deta = 0.25 * np.array([-(1-xi), -(1+xi), (1+xi), (1-xi)])
            dN_dx = (2.0 / hx) * dN_dxi
            dN_dy = (2.0 / hy) * dN_deta
            weight = hx * hy / 4.0
            Ke0 += (np.outer(dN_dx, dN_dx) + np.outer(dN_dy, dN_dy)) * weight
            fe0x += -dN_dx * weight
            fe0y += -dN_dy * weight
            dNdx_list.append(dN_dx)
            dNdy_list.append(dN_dy)
            weights.append(weight)
    return Ke0, fe0x, fe0y, np.array(dNdx_list), np.array(dNdy_list), np.array(weights)


def periodic_conn(Ny: int, Nx: int) -> np.ndarray:
    ey, ex = np.meshgrid(np.arange(Ny), np.arange(Nx), indexing="ij")
    ex = ex.ravel(); ey = ey.ravel()
    exR = (ex + 1) % Nx
    eyT = (ey + 1) % Ny
    def nid(ix, iy): return iy * Nx + ix
    conn = np.column_stack([nid(ex, ey), nid(exR, ey), nid(exR, eyT), nid(ex, eyT)]).astype(np.int64)
    return conn


def kappa_app_fem_py(kappa: np.ndarray) -> np.ndarray:
    """Periodic Q4 apparent scalar conductivity, matching kappa_app_fem.m."""
    if np.any(kappa <= 0):
        raise ValueError("kappa must be strictly positive in selected FEM window")
    Ny, Nx = kappa.shape
    Nel = Ny * Nx; Nn = Nel
    hx = 1.0 / Nx; hy = 1.0 / Ny
    conn = periodic_conn(Ny, Nx)
    kval = kappa.ravel()
    Ke0, fe0x, fe0y, dNdx, dNdy, weights = q4_scalar_local(hx, hy)
    I = np.repeat(conn, 4, axis=1).ravel()
    J = np.tile(conn, (1, 4)).ravel()
    V = (kval[:, None] * Ke0.T.reshape(1, 16)).ravel()
    K = sp.coo_matrix((V, (I, J)), shape=(Nn, Nn)).tocsr()
    K = 0.5 * (K + K.T)
    Fxvals = (kval[:, None] * fe0x[None, :]).ravel()
    Fyvals = (kval[:, None] * fe0y[None, :]).ravel()
    fx = np.bincount(conn.ravel(), weights=Fxvals, minlength=Nn)
    fy = np.bincount(conn.ravel(), weights=Fyvals, minlength=Nn)
    F = np.column_stack([fx, fy])
    area_elem = hx * hy
    m = np.bincount(conn.ravel(), weights=np.full(4*Nel, area_elem/4), minlength=Nn)
    Kaug = sp.bmat([[K, sp.csr_matrix(m[:, None])], [sp.csr_matrix(m[None, :]), sp.csr_matrix((1, 1))]], format="csr")
    Faug = np.vstack([F, np.zeros((1, 2))])
    Uaug = spla.spsolve(Kaug, Faug)
    U = np.asarray(Uaug[:Nn, :])
    # average conductivity: <kappa (e_j + grad u_j)>
    Kapp = np.zeros((2, 2), float)
    Ue_all = U[conn, :]  # Nel x 4 x 2
    for j in range(2):
        e = np.zeros(2); e[j] = 1.0
        Ue = Ue_all[:, :, j]
        for q in range(4):
            gx = Ue @ dNdx[q]
            gy = Ue @ dNdy[q]
            Kapp[0, j] += np.sum(kval * (e[0] + gx)) * weights[q]
            Kapp[1, j] += np.sum(kval * (e[1] + gy)) * weights[q]
    return Kapp

# =============================================================================
# Periodic Q4 full elasticity solver
# =============================================================================

def stiffness_tensor(E: np.ndarray, nu: np.ndarray, mode: str = "plane_strain") -> Tuple[np.ndarray, ...]:
    """Return C11,C12,C13,C21,C22,C23,C31,C32,C33 arrays in GPa."""
    if np.any(E <= 0):
        raise ValueError("E must be strictly positive")
    if mode == "plane_strain":
        if np.any(nu <= -1) or np.any(nu >= 0.5):
            raise ValueError("plane_strain requires -1 < nu < 0.5")
        c = E / ((1 + nu) * (1 - 2 * nu))
        C11 = c * (1 - nu)
        C22 = C11.copy()
        C12 = c * nu
        C33 = c * (0.5 - nu)
    elif mode == "plane_stress":
        if np.any(nu <= -1) or np.any(nu >= 1):
            raise ValueError("plane_stress requires -1 < nu < 1")
        c = E / (1 - nu**2)
        C11 = c
        C22 = c
        C12 = c * nu
        C33 = E / (2 * (1 + nu))
    else:
        raise ValueError("--elastic-mode must be plane_strain or plane_stress")
    zeros = np.zeros_like(E)
    return C11, C12, zeros, C12, C22, zeros, zeros, zeros, C33


def q4_B_matrices(hx: float, hy: float):
    g = np.array([-1.0, 1.0]) / math.sqrt(3.0)
    Bs = []
    weights = []
    for xi in g:
        for eta in g:
            dN_dxi = 0.25 * np.array([-(1-eta), (1-eta), (1+eta), -(1+eta)])
            dN_deta = 0.25 * np.array([-(1-xi), -(1+xi), (1+xi), (1-xi)])
            dN_dx = (2.0 / hx) * dN_dxi
            dN_dy = (2.0 / hy) * dN_deta
            B = np.zeros((3, 8), float)
            B[0, 0::2] = dN_dx
            B[1, 1::2] = dN_dy
            B[2, 0::2] = dN_dy
            B[2, 1::2] = dN_dx
            Bs.append(B)
            weights.append(hx * hy / 4.0)
    return Bs, np.array(weights)


def c_app_fem_py(E: np.ndarray, nu: np.ndarray, mode: str = "plane_strain") -> np.ndarray:
    """Periodic Q4 apparent 2D stiffness tensor, matching c_app_fem.m."""
    Ny, Nx = E.shape
    Nel = Ny * Nx; Nn = Nel; Ndof = 2 * Nn
    hx = 1.0 / Nx; hy = 1.0 / Ny
    conn = periodic_conn(Ny, Nx)
    edof = np.empty((Nel, 8), dtype=np.int64)
    edof[:, 0::2] = 2 * conn
    edof[:, 1::2] = 2 * conn + 1
    C11, C12, C13, C21, C22, C23, C31, C32, C33 = [x.ravel() for x in stiffness_tensor(E, nu, mode)]
    Bs, weights = q4_B_matrices(hx, hy)
    KeVals = np.zeros((Nel, 64), float)
    FeVals = np.zeros((Nel, 8, 3), float)
    for B, weight in zip(Bs, weights):
        for a in range(8):
            Ba1, Ba2, Ba3 = B[0, a], B[1, a], B[2, a]
            for b in range(8):
                Bb1, Bb2, Bb3 = B[0, b], B[1, b], B[2, b]
                idx = a * 8 + b
                KeVals[:, idx] += weight * (
                    Ba1 * (C11 * Bb1 + C12 * Bb2 + C13 * Bb3) +
                    Ba2 * (C21 * Bb1 + C22 * Bb2 + C23 * Bb3) +
                    Ba3 * (C31 * Bb1 + C32 * Bb2 + C33 * Bb3)
                )
            FeVals[:, a, 0] += -weight * (Ba1 * C11 + Ba2 * C21 + Ba3 * C31)
            FeVals[:, a, 1] += -weight * (Ba1 * C12 + Ba2 * C22 + Ba3 * C32)
            FeVals[:, a, 2] += -weight * (Ba1 * C13 + Ba2 * C23 + Ba3 * C33)
    I = np.repeat(edof, 8, axis=1).ravel()
    J = np.tile(edof, (1, 8)).ravel()
    K = sp.coo_matrix((KeVals.ravel(), (I, J)), shape=(Ndof, Ndof)).tocsr()
    K = 0.5 * (K + K.T)
    F = np.zeros((Ndof, 3), float)
    for j in range(3):
        F[:, j] = np.bincount(edof.ravel(), weights=FeVals[:, :, j].ravel(), minlength=Ndof)
    area_elem = hx * hy
    mNode = np.bincount(conn.ravel(), weights=np.full(4*Nel, area_elem/4), minlength=Nn)
    G = sp.lil_matrix((2, Ndof))
    G[0, 0::2] = mNode
    G[1, 1::2] = mNode
    G = G.tocsr()
    Kaug = sp.bmat([[K, G.T], [G, sp.csr_matrix((2, 2))]], format="csr")
    Faug = np.vstack([F, np.zeros((2, 3))])
    Daug = spla.spsolve(Kaug, Faug)
    D = np.asarray(Daug[:Ndof, :])
    # Average stiffness: Capp(:,j)=< C (e_j + B d_j) >
    Capp = np.zeros((3, 3), float)
    for j in range(3):
        De = D[edof, j]  # Nel x 8
        for B, weight in zip(Bs, weights):
            eps = De @ B.T
            eps[:, j] += 1.0
            sig1 = C11 * eps[:, 0] + C12 * eps[:, 1] + C13 * eps[:, 2]
            sig2 = C21 * eps[:, 0] + C22 * eps[:, 1] + C23 * eps[:, 2]
            sig3 = C31 * eps[:, 0] + C32 * eps[:, 1] + C33 * eps[:, 2]
            Capp[0, j] += np.sum(sig1) * weight
            Capp[1, j] += np.sum(sig2) * weight
            Capp[2, j] += np.sum(sig3) * weight
    Capp = 0.5 * (Capp + Capp.T)
    return Capp

# =============================================================================
# Statistics and plotting
# =============================================================================

def sem_or_std(vals: np.ndarray, mode: str) -> float:
    vals = np.asarray(vals, float)
    vals = vals[np.isfinite(vals)]
    if len(vals) <= 1:
        return np.nan
    s = float(np.std(vals, ddof=1))
    return s / math.sqrt(len(vals)) if mode == "sem" else s


def summarize(rows: List[Dict], Lvals: Sequence[int], mode: str) -> Dict[int, Dict[str, float]]:
    fields = ["K11", "K22", "C11", "C12", "C22", "C33", "Ex", "Ey", "nuxy", "nuyx"]
    out: Dict[int, Dict[str, float]] = {}
    for L in Lvals:
        rr = [r for r in rows if r["L"] == L]
        d: Dict[str, float] = {"n": len(rr)}
        for f in fields:
            vals = np.array([r[f] for r in rr], float)
            d[f + "_mean"] = float(np.nanmean(vals)) if len(vals) else np.nan
            d[f + "_err"] = sem_or_std(vals, mode) if len(vals) else np.nan
        out[L] = d
    return out


def save_csvs(outdir: Path, prefix: str, rows: List[Dict], summary: Dict[int, Dict[str, float]], Lvals: Sequence[int]):
    outdir.mkdir(parents=True, exist_ok=True)
    if rows:
        keys = list(rows[0].keys())
        with open(outdir / f"{prefix}_window_values.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader(); w.writerows(rows)
    keys = ["L", "L_mm", "n"]
    for f in ["K11", "K22", "C11", "C12", "C22", "C33", "Ex", "Ey", "nuxy", "nuyx"]:
        keys += [f + "_mean", f + "_err"]
    with open(outdir / f"{prefix}_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for L in Lvals:
            row = {"L": L, **summary[L]}
            row.setdefault("L_mm", np.nan)
            w.writerow(row)


def plot_err(ax, x, y, e, label, marker):
    y = np.asarray(y, float); e = np.asarray(e, float)
    low = np.minimum(e, 0.999 * np.maximum(y, EPS))
    ax.errorbar(x, y, yerr=[low, e], marker=marker, lw=1.6, ms=5, capsize=3, label=label)


def make_plots(outdir: Path, prefix: str, Lvals: Sequence[int], summary: Dict[int, Dict[str, float]],
               L_rea: int, px_mm: float, dpi: int, save_pdf: bool):
    x = np.array(Lvals, float) * px_mm
    Lrea_mm = L_rea * px_mm
    def arr(key): return np.array([summary[L][key] for L in Lvals], float)

    # Kappa
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    plot_err(ax, x, arr("K11_mean"), arr("K11_err"), r"$\kappa_x^{\rm app}$", "o")
    plot_err(ax, x, arr("K22_mean"), arr("K22_err"), r"$\kappa_y^{\rm app}$", "s")
    ax.axvline(Lrea_mm, ls="--", color="k", lw=1.4, label=fr"$L_{{\rm REA}}={L_rea}$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"Window size, $L$ [mm]")
    ax.set_ylabel(r"$\kappa_{\rm app}$ [W m$^{-1}$ K$^{-1}$]")
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=9)
    ax.text(0.04, 0.92, "(b)", transform=ax.transAxes, fontsize=15, fontweight="bold")
    fig.tight_layout(); fig.savefig(outdir / f"{prefix}_kappa_app_curves.png", dpi=dpi)
    if save_pdf: fig.savefig(outdir / f"{prefix}_kappa_app_curves.pdf")
    plt.close(fig)

    # C components
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    for key, lab, mk in [("C11", r"$C_{11}^{\rm app}$", "o"), ("C12", r"$C_{12}^{\rm app}$", "s"), ("C22", r"$C_{22}^{\rm app}$", "^"), ("C33", r"$C_{33}^{\rm app}$", "D")]:
        plot_err(ax, x, arr(key+"_mean"), arr(key+"_err"), lab, mk)
    ax.axvline(Lrea_mm, ls="--", color="k", lw=1.4)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"Window size, $L$ [mm]")
    ax.set_ylabel(r"$C_{\rm app}$ [GPa]")
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8)
    ax.text(0.04, 0.92, "(c)", transform=ax.transAxes, fontsize=15, fontweight="bold")
    fig.tight_layout(); fig.savefig(outdir / f"{prefix}_C_app_curves.png", dpi=dpi)
    if save_pdf: fig.savefig(outdir / f"{prefix}_C_app_curves.pdf")
    plt.close(fig)

    # E directional
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    plot_err(ax, x, arr("Ex_mean"), arr("Ex_err"), r"$E_x^{\rm app}$", "o")
    plot_err(ax, x, arr("Ey_mean"), arr("Ey_err"), r"$E_y^{\rm app}$", "s")
    ax.axvline(Lrea_mm, ls="--", color="k", lw=1.4)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"Window size, $L$ [mm]")
    ax.set_ylabel(r"$E_{\rm app}$ [GPa]")
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=9)
    ax.text(0.04, 0.92, "(d)", transform=ax.transAxes, fontsize=15, fontweight="bold")
    fig.tight_layout(); fig.savefig(outdir / f"{prefix}_E_app_curves.png", dpi=dpi)
    if save_pdf: fig.savefig(outdir / f"{prefix}_E_app_curves.pdf")
    plt.close(fig)

# =============================================================================
# Main
# =============================================================================

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Full kappa and elastic homogenization inside property-grid stationary masks."
    )
    p.add_argument("--kappa-maps", nargs="+", required=True, help="CSV/NPY numerical kappa maps [W m^-1 K^-1].")
    p.add_argument("--E-maps", nargs="+", required=True, help="CSV/NPY numerical Young modulus maps [GPa].")
    p.add_argument("--nu-maps", nargs="+", required=True, help="CSV/NPY numerical Poisson-ratio maps [-].")
    p.add_argument("--masks", nargs="+", default=None, help="Explicit stationary mask files, one per map.")
    p.add_argument("--mask-dir", default=None, help="Directory containing stationary masks.")
    p.add_argument("--mask-template", default="property_combined_{i}_stationary_mask.tif", help="Mask filename template using {i} and/or {stem}.")
    p.add_argument("--window-sizes", required=True, help="Comma-separated L values in property-map pixels.")
    p.add_argument("--L-rea", type=int, required=True, help="Transferred structural REA marker in property-map pixels, e.g. 204.")
    p.add_argument("--property-pixel-size-mm", type=float, default=0.01, help="Property-map pixel size in mm/pixel.")
    p.add_argument("--min-mask-fraction", type=float, default=0.98, help="Minimum fraction of LxL window inside stationary mask.")
    p.add_argument("--n-window-samples", type=int, default=4, help="Accepted windows per image and L. Full FEM is expensive; start small.")
    p.add_argument("--elastic-mode", choices=["plane_strain", "plane_stress"], default="plane_strain", help="2D elasticity assumption.")
    p.add_argument("--errorbar-mode", choices=["sem", "std"], default="sem", help="Error bars across all accepted windows/images.")
    p.add_argument("--rng-seed", type=int, default=12345, help="Random seed for window sampling.")
    p.add_argument("--output-dir", default="property_full_homogenized_property_masks", help="Output directory.")
    p.add_argument("--output-prefix", default="figure4_property_full_homogenized_property_masks", help="Output filename prefix.")
    p.add_argument("--dpi", type=int, default=600, help="PNG resolution.")
    p.add_argument("--no-pdf", action="store_true", help="Do not write PDF figures.")
    p.add_argument("--skip-elastic", action="store_true", help="Compute only kappa_app to save time.")
    return p


def main():
    args = build_argparser().parse_args()
    n = len(args.kappa_maps)
    if len(args.E_maps) != n or len(args.nu_maps) != n:
        raise ValueError("--kappa-maps, --E-maps and --nu-maps must have same length")
    Lvals = parse_int_list(args.window_sizes)
    outdir = Path(args.output_dir); outdir.mkdir(parents=True, exist_ok=True)
    stems = [Path(p).stem for p in args.kappa_maps]
    mask_paths = resolve_masks(args, n, stems)
    rng = np.random.default_rng(args.rng_seed)

    print("Full property homogenization inside property-grid stationary masks")
    print("="*88)
    print(f"Number of property sets        : {n}")
    print(f"Window sizes [pixels]         : {Lvals}")
    print(f"Property pixel size [mm/pix]  : {args.property_pixel_size_mm}")
    print(f"L_REA marker                  : {args.L_rea} px = {args.L_rea*args.property_pixel_size_mm:.4g} mm")
    print(f"Elastic mode                  : {args.elastic_mode}")
    print(f"Samples per image per L       : {args.n_window_samples}")
    print("="*88)

    rows: List[Dict] = []
    for idx in range(n):
        kappa = load_numeric_map(args.kappa_maps[idx])
        E = load_numeric_map(args.E_maps[idx])
        nu = load_numeric_map(args.nu_maps[idx])
        mask = load_mask(mask_paths[idx])
        if not (kappa.shape == E.shape == nu.shape == mask.shape):
            raise ValueError(f"Shape mismatch image {idx+1}: kappa {kappa.shape}, E {E.shape}, nu {nu.shape}, mask {mask.shape}")
        print(f"\nImage {idx+1}/{n}: shape={kappa.shape}, mask fraction={mask.mean():.3f}")
        for L in Lvals:
            wins = sample_valid_windows(mask, L, args.n_window_samples, args.min_mask_fraction, rng)
            print(f"  L={L:4d}: accepted {len(wins):3d} windows")
            for iw, (y0, x0, frac) in enumerate(wins, start=1):
                t0 = time.time()
                kw = kappa[y0:y0+L, x0:x0+L]
                Ew = E[y0:y0+L, x0:x0+L]
                nw = nu[y0:y0+L, x0:x0+L]
                Kapp = kappa_app_fem_py(kw)
                if args.skip_elastic:
                    Capp = np.full((3,3), np.nan); Ex=Ey=nuxy=nuyx=np.nan
                else:
                    Capp = c_app_fem_py(Ew, nw, args.elastic_mode)
                    try:
                        Sapp = np.linalg.inv(Capp)
                        Ex = 1.0 / Sapp[0, 0]
                        Ey = 1.0 / Sapp[1, 1]
                        nuxy = -Sapp[0, 1] / Sapp[0, 0]
                        nuyx = -Sapp[0, 1] / Sapp[1, 1]
                    except np.linalg.LinAlgError:
                        Ex=Ey=nuxy=nuyx=np.nan
                rows.append({
                    "image": idx+1, "window_index": iw, "L": L, "L_mm": L*args.property_pixel_size_mm,
                    "y0": y0, "x0": x0, "mask_fraction": frac,
                    "K11": Kapp[0,0], "K12": Kapp[0,1], "K21": Kapp[1,0], "K22": Kapp[1,1],
                    "C11": Capp[0,0], "C12": Capp[0,1], "C13": Capp[0,2],
                    "C22": Capp[1,1], "C23": Capp[1,2], "C33": Capp[2,2],
                    "Ex": Ex, "Ey": Ey, "nuxy": nuxy, "nuyx": nuyx,
                    "runtime_s": time.time()-t0,
                })
    summary = summarize(rows, Lvals, args.errorbar_mode)
    for L in Lvals:
        summary[L]["L_mm"] = L * args.property_pixel_size_mm
    save_csvs(outdir, args.output_prefix, rows, summary, Lvals)
    make_plots(outdir, args.output_prefix, Lvals, summary, args.L_rea, args.property_pixel_size_mm, args.dpi, not args.no_pdf)
    print("\nDone.")
    print(f"Outputs written to: {outdir}")

if __name__ == "__main__":
    main()
