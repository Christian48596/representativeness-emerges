#!/usr/bin/env python3
"""Publication layout for Figure 2, preserving the original numerical curves.

Usage: python regenerate_figure2_layout.py [--data-root /path/to/archive]
Data root contains BSE/1_cut.tif and reference_outputs/.
The compact numerical display data are in derived_data/figure2_plot_data.npz.
No analysis or original computational code is modified. Fixed panel boxes have
identical dimensions and the natural aspect ratio of the complete BSE image.
"""
from pathlib import Path
import argparse,json,os
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

HERE=Path(__file__).resolve().parent
REPO=HERE.parent
Image.MAX_IMAGE_PIXELS=None

def median_trend(v,width=3,log=False):
    v=np.asarray(v,dtype=float)
    x=np.where(np.isfinite(v)&(v>0),np.log10(np.maximum(v,1e-4)),np.nan) if log else v
    out=np.full(v.shape,np.nan)
    for j in range(len(v)):
        w=x[max(0,j-width//2):min(len(v),j+width//2+1)]
        good=w[np.isfinite(w)]
        if len(good):out[j]=np.median(good)
    return 10**out if log else out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--data-root',type=Path,default=Path(os.environ.get('REA_DATA_ROOT', str(REPO/'data'))))
    ap.add_argument('--data-file',type=Path,default=REPO/'derived_data/figure2_plot_data.npz')
    ap.add_argument('--output-dir',type=Path,default=REPO/'outputs/figures')
    ap.add_argument('--dpi',type=int,default=600)
    args=ap.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    d=np.load(args.data_file)
    raw=Image.open(args.data_root/'BSE/1_cut.tif').convert('L')
    W,H=raw.size
    # Display-only antialiasing to exceed the 2160-pixel panel resolution;
    # preserves the full field of view and original coordinate extent.
    raw.thumbnail((3000,3000),Image.Resampling.LANCZOS)
    image=np.asarray(raw,dtype=np.uint8)
    L=d['L'];means=d['sample_mean_Z'];eta=d['eta_C'];err=d['eta_C_sem']
    tauz=float(d['tau_Z']);tauc=float(d['tau_C']);accepted=int(d['accepted_L'])
    metric=np.where(np.isfinite(eta),np.maximum(eta,1e-4),np.nan)/tauc
    metric_err=np.where(np.isfinite(err),np.maximum(err,1e-4),np.nan)/tauc
    metric_err[L==int(d['reference_L'])]=np.nan
    smooth_mean=median_trend(means)
    smooth_metric=median_trend(metric,log=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
        'axes.labelsize':11,'axes.linewidth':.85,'xtick.labelsize':10,'ytick.labelsize':10,
        'xtick.major.width':.8,'ytick.major.width':.8,'pdf.fonttype':42,'ps.fonttype':42})
    FW,FH=10.5,8.7
    pw=3.60;ph=pw*H/W
    fig=plt.figure(figsize=(FW,FH),facecolor='white')
    rects=[[.85/FW,4.87/FH,pw/FW,ph/FH],[6.15/FW,4.87/FH,pw/FW,ph/FH],
           [.85/FW,.72/FH,pw/FW,ph/FH],[6.15/FW,.72/FH,pw/FW,ph/FH]]
    axes=[fig.add_axes(r) for r in rects]
    for ax,lab in zip(axes,'abcd'):
        ax.text(0,1.035,f'({lab})',transform=ax.transAxes,ha='left',va='bottom',
                fontsize=16,fontweight='bold',clip_on=False)
        ax.tick_params(direction='out',length=4,pad=4)
    # (a) The complete original image; square is centred exactly as in v24.
    a=axes[0]
    im=a.imshow(image,cmap='gray',origin='lower',interpolation='antialiased',extent=(-.5,W-.5,-.5,H-.5))
    x0=W//2-accepted//2;y0=H//2-accepted//2
    a.add_patch(plt.Rectangle((x0,y0),accepted,accepted,fill=False,edgecolor='white',linewidth=1.7))
    a.set_xlabel('x (pixels)');a.set_ylabel('y (pixels)')
    a.set_xticks([0,5000,10000,15000]);a.set_yticks([0,4000,8000,12000])
    cbax=fig.add_axes([4.64/FW,(4.87+.13*ph)/FH,.13/FW,.74*ph/FH])
    cb=fig.colorbar(im,cax=cbax,ticks=[0,50,100,150,200,250])
    cb.set_label(r'$Z_{\mathrm{gray}}$',fontsize=10,labelpad=5)
    cb.ax.tick_params(labelsize=9,pad=3,length=3)
    # (b) Original finite-window means and three-point running median.
    b=axes[1]
    refmean=means[L==accepted][0]
    b.axhspan(refmean*(1-tauz),refmean*(1+tauz),color='C0',alpha=.08,linewidth=0,label=r'$\pm\tau_Z$ band')
    b.plot(L,means,'o-',color='C0',alpha=.35,lw=1.15,ms=4,label=r'$\overline{Z}_L$')
    b.plot(L,smooth_mean,'o-',color='C0',lw=1.75,ms=4.4,label='Median trend')
    b.axvline(accepted,color='black',ls='--',lw=1.5,alpha=.9,label=rf'$L_{{\mathrm{{REA}}}}={accepted}$')
    b.set_xlabel(r'Window size $L$ (pixels)');b.set_ylabel(r'Ensemble mean $\overline{Z}_L$')
    b.grid(True,alpha=.20)
    b.legend(frameon=False,fontsize=8.5,loc='lower left',handlelength=2.2)
    # (c) Exactly the original raw ratios, median guide and post-selection SEM.
    c=axes[2]
    positive=metric[np.isfinite(metric)&(metric>0)]
    ymin=max(min(positive.min()/2,.8),1e-4);ymax=max(positive.max()*1.8,1.35)
    c.fill_between([L.min(),L.max()],ymin,1,color='.92',alpha=.70,zorder=0)
    c.plot(L,metric,'s-',color='C1',alpha=.35,lw=1.05,ms=4,label=r'$\eta_{\widehat C}/\tau_{\widehat C}$')
    c.plot(L,smooth_metric,'s-',color='C1',lw=1.8,ms=4.6,label='Median trend')
    take=(L>=accepted)&np.isfinite(metric)&np.isfinite(metric_err)
    c.errorbar(L[take],metric[take],yerr=metric_err[take],fmt='none',ecolor='C1',elinewidth=1,
               capsize=3,capthick=1,zorder=2.5,alpha=.9)
    c.set_yscale('log');c.axhline(1,color='black',ls=':',lw=1.1)
    c.axvline(accepted,color='black',ls='--',lw=1.5,alpha=.9,label=rf'$L_{{\mathrm{{REA}}}}={accepted}$')
    c.set_ylim(ymin,ymax);c.set_xlim(L.min()-.02*np.ptp(L),L.max()+.07*np.ptp(L))
    c.set_xlabel(r'Window size $L$ (pixels)');c.set_ylabel('Structural residual / tolerance')
    c.grid(True,which='major',alpha=.24);c.grid(True,which='minor',alpha=.10)
    c.legend(frameon=True,facecolor='white',edgecolor='none',framealpha=.94,fontsize=8.5,
             loc='upper right',handlelength=2.2)
    for ax in [b,c]:ax.set_xticks([256,512,1024,1536,2048])
    # (d) Original ensemble spectra and SEM shading, each peak normalized.
    a4=axes[3];xs=[];ys=[]
    for n in [256,512,1024,2048]:
        k=d[f'k_{n}'];p=d[f'p_{n}'];perr=d[f'sem_{n}']
        good=(k>0)&(p>0)&np.isfinite(k)&np.isfinite(p)
        norm=np.nanmax(p[good])+1e-30;y=p[good]/norm;ye=perr[good]/norm
        lo=np.maximum(y-ye,1e-30);hi=y+ye
        a4.loglog(k[good],y,lw=1.4,label=rf'$L={n}$')
        a4.fill_between(k[good],lo,hi,alpha=.10,lw=0)
        xs.append(k[good]);ys.extend([y,lo,hi])
    xx=np.concatenate(xs);yy=np.concatenate(ys)
    a4.set_xlim(xx.min()/1.15,xx.max()*1.15)
    a4.set_ylim(max(yy[yy>0].min()/1.8,1e-8),min(yy.max()*1.45,2.0))
    a4.set_xlabel(r'Wavenumber $k$ (pixel$^{-1}$)')
    a4.set_ylabel(r'Normalized spectrum $\widehat C_L(k)/\max_k\widehat C_L(k)$')
    a4.grid(True,which='both',alpha=.20)
    a4.legend(frameon=False,fontsize=8.3,ncol=2,loc='upper right',handlelength=1.8,columnspacing=1.1)
    # Verify actual rendered axes sizes, including image aspect constraints.
    fig.canvas.draw()
    boxes=np.array([[ax.get_window_extent().width,ax.get_window_extent().height] for ax in axes])
    assert np.max(np.abs(boxes-boxes[0]))<1e-7,boxes
    assert all(not ax.get_title() for ax in axes)
    out=args.output_dir/'Figure_2_revised'
    fig.savefig(out.with_suffix('.png'),dpi=args.dpi,facecolor='white')
    fig.savefig(out.with_suffix('.pdf'),dpi=args.dpi,facecolor='white')
    fig.savefig(args.output_dir/'Figure_2_preview.png',dpi=145,facecolor='white')
    audit={'panel_width_inches':pw,'panel_height_inches':ph,'identical_panel_boxes':True,
           'original_image_shape':[H,W],'display_image_shape':list(image.shape),'image_source':str(args.data_root/'BSE/1_cut.tif'),
           'accepted_L':accepted,'reference_L':int(d['reference_L']),
           'raw_sample_mean_Z':means.tolist(),'raw_eta_C':eta.tolist(),'eta_C_sem':err.tolist(),
           'display_only_changes':['Equal-sized fixed panel boxes','Bold external upper-left labels',
                'Font/legend placement','BSE raster antialiased for display at 3000px; complete field of view retained','Corrected spectral axis labels for existing log axes; no numerical transform'],
           'numerical_values_unchanged':True}
    (args.output_dir/'Figure_2_layout_audit.json').write_text(json.dumps(audit,indent=2))
    plt.close(fig)
    print(f'Saved {out}.png and .pdf; all panel boxes {pw:g} x {ph:g} inches')

if __name__=='__main__':main()
