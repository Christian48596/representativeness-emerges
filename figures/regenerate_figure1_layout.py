#!/usr/bin/env python3
"""Plot an image, admissibility mask and detector diagnostics as equal panels.

No statistical calculation is performed. Explicit input paths support arbitrary
image names; without them, the bundled study's image-4 paths are used.
"""
from pathlib import Path
import argparse, csv, json, os
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FormatStrFormatter

REPO = Path(__file__).resolve().parents[1]
Image.MAX_IMAGE_PIXELS = None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, default=Path(os.environ.get('REA_DATA_ROOT', str(REPO/'data'))))
    parser.add_argument('--image-file', type=Path, help='Representative image; any image filename is accepted.')
    parser.add_argument('--mask-file', type=Path, help='Pixel mask matching the image dimensions.')
    parser.add_argument('--windows-file', type=Path, help='Detector window-diagnostic CSV for this image.')
    parser.add_argument('--output-dir', type=Path, default=REPO/'outputs/figures')
    parser.add_argument('--dpi', type=int, default=600)
    args = parser.parse_args()
    if args.dpi < 1:
        parser.error('--dpi must be positive')
    src = args.data_root/'reference_outputs/stationary_domains_FULLRES_R2048_STRICT'
    image_path = args.image_file or args.data_root/'BSE/4_cut.tif'
    mask_path = args.mask_file or src/f'{image_path.stem}_stationary_mask.tif'
    windows_path = args.windows_file or src/f'{image_path.stem}_stationary_windows.csv'
    for path in (image_path, mask_path, windows_path):
        if not path.is_file():
            parser.error(f'Input file not found: {path}. Supply --image-file, --mask-file and --windows-file for custom data.')
    with Image.open(image_path) as im:
        width, height = im.size
        display = im.convert('L')
        display.thumbnail((3000, 3000), Image.Resampling.LANCZOS)
        field = np.asarray(display)
    with Image.open(mask_path) as im:
        if im.size != (width, height):
            parser.error('The image and mask must have identical pixel dimensions.')
        raw_mask = np.asarray(im.convert('L'))
        fraction = float(np.count_nonzero(raw_mask)/raw_mask.size)
        mask = np.asarray(im.convert('L').resize(display.size, Image.Resampling.NEAREST)) > 0
    with windows_path.open(newline='') as handle:
        reader = csv.DictReader(handle)
        required = {'iy', 'ix', 'mean', 'mu_residual', 'std_residual', 'stationary_clean'}
        if not required.issubset(reader.fieldnames or []):
            parser.error('Window CSV must contain: '+', '.join(sorted(required)))
        rows = list(reader)
    if not rows:
        parser.error('The window CSV is empty; no screening grid can be displayed.')
    ny = max(int(row['iy']) for row in rows)+1
    nx = max(int(row['ix']) for row in rows)+1
    if min(ny, nx) < 1:
        parser.error('Window-grid indices must be nonnegative.')
    means = np.full((ny, nx), np.nan)
    residual = np.full_like(means, np.nan)
    clean = np.zeros((ny, nx), dtype=bool)
    seen = set()
    for row in rows:
        iy, ix = int(row['iy']), int(row['ix'])
        if iy < 0 or ix < 0 or (iy, ix) in seen:
            parser.error('Window-grid indices must be unique and nonnegative.')
        seen.add((iy, ix))
        means[iy, ix] = float(row['mean'])
        residual[iy, ix] = max(float(row['mu_residual']), float(row['std_residual']))
        clean[iy, ix] = int(row['stationary_clean']) > 0
    if len(seen) != ny*nx:
        parser.error('The window CSV does not contain a complete rectangular grid.')
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':11, 'axes.labelsize':12,
        'xtick.labelsize':10.5, 'ytick.labelsize':10.5, 'axes.linewidth':.8,
        'pdf.fonttype':42, 'ps.fonttype':42, 'savefig.facecolor':'white'})
    fw, pw = 10.5, 3.6
    ph = pw*height/width
    bottom, row_gap, top_margin = .74, 1.02, .65
    upper = bottom+ph+row_gap
    fh = upper+ph+top_margin
    fig = plt.figure(figsize=(fw, fh), facecolor='white')
    axes = [fig.add_axes([x/fw,y/fh,pw/fw,ph/fh]) for x,y in
            [(.82,upper),(5.91,upper),(.82,bottom),(5.91,bottom)]]
    extent = (-.5,width-.5,height-.5,-.5)
    image = axes[0].imshow(field,cmap='viridis',vmin=0,vmax=255,origin='upper',extent=extent,aspect='auto',interpolation='nearest',rasterized=True)
    axes[1].imshow(field,cmap='viridis',vmin=0,vmax=255,origin='upper',extent=extent,aspect='auto',interpolation='nearest',rasterized=True)
    overlay = np.zeros((*mask.shape,4),dtype=np.float32)
    overlay[...,0] = 1.
    overlay[...,3] = mask*.35
    axes[1].imshow(overlay,origin='upper',extent=extent,aspect='auto',interpolation='nearest',rasterized=True)
    for ax in axes[:2]:
        ax.set_xlim(-.5,width-.5); ax.set_ylim(height-.5,-.5)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4,integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4,integer=True))
        ax.set_xlabel('x (pixels)'); ax.set_ylabel('y (pixels)')
    diagnostics = []
    for ax, values in zip(axes[2:], [means,residual]):
        im = ax.imshow(values,cmap='viridis',origin='upper',aspect='auto',interpolation='nearest',rasterized=True)
        if min(clean.shape) >= 2 and clean.any() and not clean.all():
            ax.contour(clean.astype(float),levels=[.5],colors='white',linewidths=.9)
        ax.set_xlim(-.5,nx-.5); ax.set_ylim(ny-.5,-.5)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5,integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5,integer=True))
        ax.set_xlabel('Window-grid x'); ax.set_ylabel('Window-grid y')
        diagnostics.append(im)
    for ax,im in [(axes[0],image),(axes[2],diagnostics[0]),(axes[3],diagnostics[1])]:
        box = ax.get_position()
        cax = fig.add_axes([box.x1+.11/fw,box.y0+.06*box.height,.13/fw,.88*box.height])
        cb = fig.colorbar(im,cax=cax)
        cb.ax.tick_params(labelsize=9.5,pad=3,width=.7,length=3)
        cb.outline.set_linewidth(.7)
        if ax is axes[3]: cb.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    for tag,ax in zip('abcd',axes):
        ax.text(0,1.035,f'({tag})',transform=ax.transAxes,fontsize=16,fontweight='bold',ha='left',va='bottom',clip_on=False)
        ax.tick_params(direction='out',length=4,width=.8,pad=4)
    fig.canvas.draw()
    dimensions = np.array([ax.get_position().bounds for ax in axes])*[fw,fh,fw,fh]
    assert np.allclose(dimensions[:,2:],[pw,ph],atol=1e-10)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    for suffix in ('png','pdf'):
        fig.savefig(args.output_dir/f'Figure_1_revised.{suffix}',dpi=args.dpi)
    fig.savefig(args.output_dir/'Figure_1_preview.png',dpi=min(150,args.dpi))
    audit = {'image_shape':[height,width],'mask_fraction':fraction,'coarse_grid':[ny,nx],
        'panel_size_inches':[pw,ph],'all_panel_dimensions_equal':True,
        'image_source':str(image_path),'mask_source':str(mask_path),'windows_source':str(windows_path),
        'numerical_values_unchanged':True}
    (args.output_dir/'Figure_1_layout_audit.json').write_text(json.dumps(audit,indent=2,allow_nan=False)+'\n')
    plt.close(fig)
    print(json.dumps(audit,indent=2))


if __name__ == '__main__':
    main()
