#!/usr/bin/env python3
"""Replot Figure 1 from archived raw BSE image, mask and exact diagnostic CSV.
Only presentation changes. No REA/stationarity calculation is performed.
Usage: python regenerate_figure1_layout.py --data-root /path/to/archive
"""
from pathlib import Path
import argparse,csv,json,os
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import FormatStrFormatter

REPO=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser()
p.add_argument('--data-root',type=Path,default=Path(os.environ.get('REA_DATA_ROOT', str(REPO/'data'))))
p.add_argument('--output-dir',type=Path,default=REPO/'outputs/figures')
a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
root=a.data_root
src=root/'reference_outputs/stationary_domains_FULLRES_R2048_STRICT'
Image.MAX_IMAGE_PIXELS=None
with Image.open(root/'BSE/4_cut.tif') as im:
    W,H=im.size
    # Downsample only the display image, to above the final output panel resolution.
    display=im.convert('L');display.thumbnail((3000,3000),Image.Resampling.LANCZOS)
    field=np.asarray(display)
with Image.open(src/'4_cut_stationary_mask.tif') as im:
    assert im.size==(W,H)
    rawmask=np.asarray(im)
    fraction=np.count_nonzero(rawmask)/rawmask.size
    # Nearest-neighbour rendering preserves the discrete mask assignment.
    mask=np.asarray(im.resize(display.size,Image.Resampling.NEAREST))>0
with (src/'4_cut_stationary_windows.csv').open() as f:rows=list(csv.DictReader(f))
ny=max(int(r['iy']) for r in rows)+1;nx=max(int(r['ix']) for r in rows)+1
means=np.empty((ny,nx));residual=np.empty_like(means);clean=np.empty((ny,nx),dtype=bool)
for r in rows:
    iy,ix=int(r['iy']),int(r['ix'])
    means[iy,ix]=float(r['mean'])
    residual[iy,ix]=max(float(r['mu_residual']),float(r['std_residual']))
    clean[iy,ix]=int(r['stationary_clean'])>0
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':12,
    'xtick.labelsize':10.5,'ytick.labelsize':10.5,'axes.linewidth':.8,
    'pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'white'})
fw,fh=10.5,8.7
pw=3.6;ph=pw*H/W
fig=plt.figure(figsize=(fw,fh),facecolor='white')
axes=[]
for x,y in [(0.82,4.95),(5.91,4.95),(0.82,.74),(5.91,.74)]:
    axes.append(fig.add_axes([x/fw,y/fh,pw/fw,ph/fh]))
# Image-pixel aspect is preserved exactly by matching the physical axes ratio.
extent=(-.5,W-.5,H-.5,-.5)
im0=axes[0].imshow(field,cmap='viridis',vmin=0,vmax=255,origin='upper',extent=extent,aspect='auto',interpolation='nearest',rasterized=True)
axes[1].imshow(field,cmap='viridis',vmin=0,vmax=255,origin='upper',extent=extent,aspect='auto',interpolation='nearest',rasterized=True)
overlay=np.zeros((*mask.shape,4),dtype=np.float32);overlay[...,0]=1.;overlay[...,3]=mask*.35
axes[1].imshow(overlay,origin='upper',extent=extent,aspect='auto',interpolation='nearest',rasterized=True)
for ax in axes[:2]:
    ax.set_xlim(-.5,W-.5);ax.set_ylim(H-.5,-.5)
    ax.set_xticks([0,5000,10000,15000]);ax.set_yticks([0,4000,8000,12000])
    ax.set_xlabel('x (pixels)');ax.set_ylabel('y (pixels)')
images=[]
for ax,arr in zip(axes[2:],[means,residual]):
    im=ax.imshow(arr,cmap='viridis',origin='upper',aspect='auto',interpolation='nearest',rasterized=True)
    ax.contour(clean.astype(float),levels=[.5],colors='white',linewidths=.9)
    ax.set_xlim(-.5,nx-.5);ax.set_ylim(ny-.5,-.5)
    ax.set_xticks([0,5,10,15,20,25]);ax.set_yticks([0,5,10,15,20])
    ax.set_xlabel('Window-grid x');ax.set_ylabel('Window-grid y')
    images.append(im)
# Colorbars have their own axes and cannot change panel dimensions.
for ax,im,label in [(axes[0],im0,'Gray level'),(axes[2],images[0],'Local mean'),(axes[3],images[1],'Normalized residual')]:
    box=ax.get_position();cax=fig.add_axes([box.x1+.11/fw,box.y0+.06*box.height,.13/fw,.88*box.height])
    cb=fig.colorbar(im,cax=cax);cb.ax.tick_params(labelsize=9.5,pad=3,width=.7,length=3)
    cb.outline.set_linewidth(.7)
    if ax is axes[3]:cb.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
for tag,ax in zip('abcd',axes):
    ax.text(0,1.035,f'({tag})',transform=ax.transAxes,fontsize=16,fontweight='bold',ha='left',va='bottom',clip_on=False)
    ax.tick_params(direction='out',length=4,width=.8,pad=4)
fig.canvas.draw()
dims=[np.array(ax.get_position().bounds)*[fw,fh,fw,fh] for ax in axes]
assert np.allclose(np.array(dims)[:,2:],[pw,ph],atol=1e-10)
assert all(not ax.get_title() for ax in axes)
for suffix in ['png','pdf']:
    fig.savefig(a.output_dir/f'Figure_1_revised.{suffix}',dpi=600)
fig.savefig(a.output_dir/'Figure_1_preview.png',dpi=150)
summary={'image_shape':[H,W],'mask_fraction':fraction,'coarse_grid':[ny,nx],
    'panel_size_inches':[pw,ph],'all_panel_dimensions_equal':True,
    'data_source':'Archived BSE/4_cut.tif, 4_cut_stationary_mask.tif, 4_cut_stationary_windows.csv',
    'changes':'Presentation only: no titles, bold external labels, equal panel axes, separate colorbars.'}
(a.output_dir/'Figure_1_layout_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
plt.close(fig)
