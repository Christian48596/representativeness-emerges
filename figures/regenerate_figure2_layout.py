#!/usr/bin/env python3
"""Plot structural convergence for any ensemble size and support grid.

The NPZ contains numerical ensemble curves, numeric accepted_L (NaN for no
selection), reference_L, and optional k_L/p_L/sem_L spectra. An explicit image
path selects the representative field. No statistical calculation is performed.
"""
from pathlib import Path
import argparse, json, os, re
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

REPO=Path(__file__).resolve().parents[1]
Image.MAX_IMAGE_PIXELS=None


def median_trend(values,width=3,log=False):
    values=np.asarray(values,dtype=float)
    transformed=np.where(np.isfinite(values)&(values>0),np.log10(np.maximum(values,1e-4)),np.nan) if log else values
    output=np.full(values.shape,np.nan)
    for j in range(len(values)):
        window=transformed[max(0,j-width//2):min(len(values),j+width//2+1)]
        good=window[np.isfinite(window)]
        if len(good): output[j]=np.median(good)
    return 10**output if log else output


def finite_json(value):
    if isinstance(value,dict): return {k:finite_json(v) for k,v in value.items()}
    if isinstance(value,list): return [finite_json(v) for v in value]
    if isinstance(value,float) and not np.isfinite(value): return None
    return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path,default=Path(os.environ.get('REA_DATA_ROOT',str(REPO/'data'))))
    parser.add_argument('--image-file',type=Path,help='Representative image for the first panel, with arbitrary filename.')
    parser.add_argument('--data-file',type=Path,default=REPO/'derived_data/figure2_plot_data.npz')
    parser.add_argument('--output-dir',type=Path,default=REPO/'outputs/figures')
    parser.add_argument('--dpi',type=int,default=600)
    parser.add_argument('--pixel-size-mm',type=float,help='Optional physical calibration; NPZ wavenumbers must be in inverse pixels.')
    args=parser.parse_args()
    if args.dpi<1: parser.error('--dpi must be positive')
    if args.pixel_size_mm is not None and (not np.isfinite(args.pixel_size_mm) or args.pixel_size_mm<=0):
        parser.error('--pixel-size-mm must be positive and finite')
    if not args.data_file.is_file(): parser.error(f'Numerical input not found: {args.data_file}')
    image_path=args.image_file or args.data_root/'BSE/1_cut.tif'
    if not image_path.is_file(): parser.error(f'Image not found: {image_path}; use --image-file for custom data.')
    with np.load(args.data_file,allow_pickle=False) as src:
        required={'L','sample_mean_Z','eta_C','eta_C_sem','tau_Z','tau_C','accepted_L','reference_L'}
        if not required.issubset(src.files): parser.error('NPZ lacks fields: '+', '.join(sorted(required-set(src.files))))
        data={key:src[key] for key in src.files}
    supports=np.asarray(data['L'],dtype=float)
    means=np.asarray(data['sample_mean_Z'],dtype=float)
    residual=np.asarray(data['eta_C'],dtype=float)
    errors=np.asarray(data['eta_C_sem'],dtype=float)
    if supports.ndim!=1 or len(supports)<1 or any(v.shape!=supports.shape for v in (means,residual,errors)):
        parser.error('L and the three numerical curves must be matching nonempty 1D arrays.')
    if not np.all(np.isfinite(supports)&(supports>0)) or not np.all(np.diff(supports)>0):
        parser.error('Support sizes must be finite, positive and strictly increasing.')
    if np.any(errors[np.isfinite(errors)] < 0) or np.any(residual[np.isfinite(residual)] < 0):
        parser.error('Spectral residuals and standard errors must be nonnegative or undefined.')
    tauz=float(data['tau_Z']); tauc=float(data['tau_C'])
    if not np.isfinite(tauz) or tauz<=0 or not np.isfinite(tauc) or tauc<=0:
        parser.error('The mean and spectral tolerances must be positive and finite.')
    accepted_value=float(data['accepted_L']); reference=float(data['reference_L'])
    accepted=accepted_value if np.isfinite(accepted_value) and accepted_value>0 else None
    if accepted is not None and not np.any(supports==accepted):
        parser.error('accepted_L must be a tested support, or numeric NaN for no selection.')
    if not np.isfinite(reference) or reference<=0: parser.error('reference_L must be positive and finite.')
    calibration=args.pixel_size_mm
    if calibration is None and 'pixel_size_mm' in data:
        candidate=float(data['pixel_size_mm'])
        if np.isfinite(candidate) and candidate>0: calibration=candidate
    scale=calibration or 1.
    length_unit='mm' if calibration else 'pixels'
    x=supports*scale
    with Image.open(image_path) as source:
        raw=source.convert('L'); width,height=raw.size
        raw.thumbnail((3000,3000),Image.Resampling.LANCZOS)
        image=np.asarray(raw,dtype=np.uint8)
    if accepted is not None and accepted>min(width,height):
        parser.error('The selected support is larger than the representative image; choose an image from the same ensemble.')
    metric=np.where(np.isfinite(residual),np.maximum(residual,1e-4),np.nan)/tauc
    metric_errors=np.where(np.isfinite(errors),np.maximum(errors,1e-4),np.nan)/tauc
    metric_errors[supports==reference]=np.nan
    smooth_mean=median_trend(means)
    smooth_metric=median_trend(metric,log=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':11,
        'axes.linewidth':.85,'xtick.labelsize':10,'ytick.labelsize':10,
        'xtick.major.width':.8,'ytick.major.width':.8,'pdf.fonttype':42,'ps.fonttype':42})
    fw,pw=10.8,3.6
    ph=pw*height/width
    bottom,row_gap,top_margin=.72,1.0,.65
    upper=bottom+ph+row_gap; fh=upper+ph+top_margin
    fig=plt.figure(figsize=(fw,fh),facecolor='white')
    axes=[fig.add_axes([left/fw,y/fh,pw/fw,ph/fh]) for left,y in
          [(1.15,upper),(6.45,upper),(1.15,bottom),(6.45,bottom)]]
    for ax,tag in zip(axes,'abcd'):
        ax.text(0,1.035,f'({tag})',transform=ax.transAxes,ha='left',va='bottom',fontsize=16,fontweight='bold',clip_on=False)
        ax.tick_params(direction='out',length=4,pad=4)
    a,b,c,d=axes
    im=a.imshow(image,cmap='gray',origin='lower',interpolation='antialiased',extent=(-.5,width-.5,-.5,height-.5),aspect='auto')
    if accepted is not None:
        x0=width//2-int(accepted)//2; y0=height//2-int(accepted)//2
        a.add_patch(plt.Rectangle((x0,y0),accepted,accepted,fill=False,edgecolor='white',linewidth=1.7))
    a.set_xlabel('x (pixels)'); a.set_ylabel('y (pixels)')
    a.xaxis.set_major_locator(MaxNLocator(nbins=4,integer=True)); a.yaxis.set_major_locator(MaxNLocator(nbins=4,integer=True))
    cbax=fig.add_axes([4.94/fw,(upper+.13*ph)/fh,.13/fw,.74*ph/fh])
    cb=fig.colorbar(im,cax=cbax)
    cb.set_label(r'$Z_{\mathrm{gray}}$',fontsize=10,labelpad=5)
    cb.ax.tick_params(labelsize=9,pad=3,length=3)
    band_support=accepted if accepted is not None else reference
    band_values=means[supports==band_support]
    if band_values.size and np.isfinite(band_values[0]):
        center=band_values[0]; limits=sorted((center*(1-tauz),center*(1+tauz)))
        b.axhspan(*limits,color='C0',alpha=.08,linewidth=0,label=r'$\pm\tau_Z$ display band')
    b.plot(x,means,'o-',color='C0',alpha=.35,lw=1.15,ms=4,label=r'$\overline{Z}_L$')
    b.plot(x,smooth_mean,'o-',color='C0',lw=1.75,ms=4.4,label='Median trend')
    if accepted is not None:
        b.axvline(accepted*scale,color='black',ls='--',lw=1.5,alpha=.9,label=rf'$L_{{\mathrm{{REA}}}}={accepted*scale:g}$ {length_unit}')
    else:
        status='No support selected' if np.isfinite(means).any() else 'No admissible common supports'
        b.text(.03,.97,status,transform=b.transAxes,ha='left',va='top',fontsize=9)
    b.set_xlabel(f'Window side length L ({length_unit})'); b.set_ylabel(r'Ensemble mean $\overline{Z}_L$')
    b.grid(True,alpha=.20); b.legend(frameon=False,fontsize=8.1,loc='lower left',handlelength=2.2)
    positive=metric[np.isfinite(metric)&(metric>0)]
    ymin=max(min(float(positive.min())/2,.8),1e-4) if positive.size else .1
    ymax=max(float(positive.max())*1.8,1.35) if positive.size else 10.
    c.fill_between([x.min(),x.max()],ymin,1,color='.92',alpha=.70,zorder=0)
    c.plot(x,metric,'s-',color='C1',alpha=.35,lw=1.05,ms=4,label=r'$\eta_{\widehat C}/\tau_{\widehat C}$')
    c.plot(x,smooth_metric,'s-',color='C1',lw=1.8,ms=4.6,label='Median trend')
    take=np.isfinite(metric)&np.isfinite(metric_errors)
    if accepted is not None: take &= supports>=accepted
    c.errorbar(x[take],metric[take],yerr=metric_errors[take],fmt='none',ecolor='C1',elinewidth=1,capsize=3,capthick=1,zorder=2.5,alpha=.9)
    c.set_yscale('log'); c.axhline(1,color='black',ls=':',lw=1.1)
    if accepted is not None:
        c.axvline(accepted*scale,color='black',ls='--',lw=1.5,alpha=.9,label=rf'$L_{{\mathrm{{REA}}}}={accepted*scale:g}$ {length_unit}')
    c.set_ylim(ymin,ymax)
    span=float(np.ptp(x)) or max(float(x[0]),1.)
    c.set_xlim(x.min()-.02*span,x.max()+.07*span)
    c.set_xlabel(f'Window side length L ({length_unit})'); c.set_ylabel('Structural residual / tolerance')
    c.grid(True,which='major',alpha=.24); c.grid(True,which='minor',alpha=.10)
    c.legend(frameon=True,facecolor='white',edgecolor='none',framealpha=.94,fontsize=8.1,loc='upper right',handlelength=2.2)
    tick_indices=np.unique(np.linspace(0,len(x)-1,min(5,len(x))).round().astype(int))
    for ax in (b,c): ax.set_xticks(x[tick_indices])
    spectra=sorted((int(match.group(1)),key) for key in data for match in [re.fullmatch(r'k_(\d+)',key)] if match)
    xs=[]; ys=[]; displayed=[]
    for n,key in spectra:
        if f'p_{n}' not in data: parser.error(f'Spectrum {key} is missing p_{n}.')
        k=np.asarray(data[key],dtype=float); power=np.asarray(data[f'p_{n}'],dtype=float)
        sem=np.asarray(data.get(f'sem_{n}',np.zeros_like(power)),dtype=float)
        if k.ndim!=1 or k.shape!=power.shape or sem.shape!=power.shape:
            parser.error(f'Spectrum arrays for L={n} must be matching 1D arrays.')
        good=(k>0)&(power>0)&np.isfinite(k)&np.isfinite(power)
        if not good.any(): continue
        normalization=float(np.max(power[good]))+1e-30
        xx=k[good]/scale; yy=power[good]/normalization
        uncertainty=np.where(np.isfinite(sem[good]),np.maximum(sem[good],0.),0.)/normalization
        lower=np.maximum(yy-uncertainty,1e-30); higher=yy+uncertainty
        d.loglog(xx,yy,lw=1.4,label=rf'$L={n*scale:g}$')
        d.fill_between(xx,lower,higher,alpha=.10,lw=0)
        xs.append(xx); ys.extend([yy,lower,higher]); displayed.append(n)
    if xs:
        xx=np.concatenate(xs); yy=np.concatenate(ys)
        d.set_xlim(xx.min()/1.15,xx.max()*1.15)
        d.set_ylim(max(yy[yy>0].min()/1.8,1e-8),max(1.5,min(yy.max()*1.45,2.0)))
        d.legend(frameon=False,fontsize=8.3,ncol=2,loc='upper right',handlelength=1.8,columnspacing=1.1)
    else:
        d.text(.5,.5,'No positive finite spectra available',transform=d.transAxes,ha='center',va='center',fontsize=9,wrap=True)
    d.set_xlabel(r'Wavenumber $k$ (mm$^{-1}$)' if calibration else r'Wavenumber $k$ (pixel$^{-1}$)')
    d.set_ylabel(r'Normalized spectrum $\widehat C_L(k)/\max_k\widehat C_L(k)$')
    d.grid(True,which='both',alpha=.20)
    fig.canvas.draw()
    boxes=np.array([[ax.get_window_extent().width,ax.get_window_extent().height] for ax in axes])
    assert np.max(np.abs(boxes-boxes[0]))<1e-7,boxes
    args.output_dir.mkdir(parents=True,exist_ok=True)
    out=args.output_dir/'Figure_2_revised'
    fig.savefig(out.with_suffix('.png'),dpi=args.dpi,facecolor='white')
    fig.savefig(out.with_suffix('.pdf'),dpi=args.dpi,facecolor='white')
    fig.savefig(args.output_dir/'Figure_2_preview.png',dpi=min(145,args.dpi),facecolor='white')
    audit={'panel_width_inches':pw,'panel_height_inches':ph,'identical_panel_boxes':True,
        'original_image_shape':[height,width],'display_image_shape':list(image.shape),'image_source':str(image_path),
        'accepted_L':accepted,'reference_L':reference,'pixel_size_mm':calibration,
        'raw_sample_mean_Z':means.tolist(),'raw_eta_C':residual.tolist(),'eta_C_sem':errors.tolist(),
        'spectrum_supports':displayed,'numerical_values_unchanged':True}
    if 'n_images' in data: audit['n_images']=int(data['n_images'])
    (args.output_dir/'Figure_2_layout_audit.json').write_text(json.dumps(finite_json(audit),indent=2,allow_nan=False)+'\n')
    plt.close(fig)
    print(f'Saved {out}.png and .pdf; selected support: {accepted}; equal panel boxes {pw:g} x {ph:g} inches')


if __name__=='__main__':
    main()
