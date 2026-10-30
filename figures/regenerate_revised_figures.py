#!/usr/bin/env python3
"""Regenerate revised panels from the archived numerical inputs and cell outputs.

No homogenization or REA code is changed or rerun. Figure 3(a) pools pixels of
all seven complete QEMSCAN maps (area weighting). Figure 3(d) inverts each
archived full symmetric plane-strain stiffness matrix before averaging the
resulting moduli; pooled-window SEM reproduces the original display convention.

Usage: python regenerate_revised_figures.py
To recompute counts from original maps: --data-root /path/to/extracted/archive
The default uses the compact verified inputs bundled in derived_data/.
"""
from pathlib import Path
import argparse, csv, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser()
ap.add_argument('--data-root','--archive',dest='archive',type=Path,default=Path(os.environ['REA_DATA_ROOT']) if os.environ.get('REA_DATA_ROOT') else None,help='Optional extracted archive (or REA_DATA_ROOT); otherwise bundled verified inputs.')
ap.add_argument('--output-dir','--output',dest='output',type=Path,default=REPO/'outputs/revised_figures')
a=ap.parse_args()
out=a.output
bundled=REPO/'derived_data'
(out/'IMG').mkdir(parents=True,exist_ok=True)
(out/'derived_data').mkdir(parents=True,exist_ok=True)
# Class identities are jointly determined by all three archived property fields.
classes=[('Quartz',6.5,94,.08),('K-feldspar',2.3,70,.25),('Plagioclase',2.,76,.26),
 ('Mica',1.,60,.25),('Calcite',3.6,72,.31),('Dolomite',5.,94,.28),
 ('Clay minerals',1.5,20,.3),('Gypsum/anhydrite',2.5,45,.3),
 ('Heavy minerals',5.,100,.25),('Pores/background',.026,.01,.2),
 ('Other/unknown',None,.001,.2)]
count_rows=[]; total_counts=np.zeros(len(classes),dtype=np.int64)
if a.archive is not None:
 for i in range(1,8):
  k,E,nu=[np.loadtxt(a.archive/'QEMSCAN'/f'{name}Map_{i}.csv',delimiter=',') for name in ('kappa','E','nu')]
  assert k.shape==E.shape==nu.shape
  seen=np.zeros(k.shape,dtype=np.uint8)
  for j,(name,k0,e0,n0) in enumerate(classes):
   mask=np.isclose(E,e0,rtol=0,atol=1e-8)&np.isclose(nu,n0,rtol=0,atol=1e-8)
   if k0 is not None: mask &= np.isclose(k,k0,rtol=0,atol=1e-8)
   count=int(mask.sum()); seen += mask
   total_counts[j]+=count
   count_rows.append({'image':i,'class':name,'pixel_count':count,'total_pixels':k.size,'area_fraction':count/k.size})
  assert np.all(seen==1), f'Unclassified/ambiguous pixels in image {i}'
else:
 with (bundled/'phase_area_fractions.csv').open() as f:
  for r in csv.DictReader(f):
   count_rows.append({'image':int(r['image']),'class':r['class'],'pixel_count':int(r['pixel_count']),'total_pixels':int(r['total_pixels']),'area_fraction':float(r['area_fraction'])})
 for j,(name,*_) in enumerate(classes):total_counts[j]=sum(r['pixel_count'] for r in count_rows if r['class']==name)
with (out/'derived_data'/'phase_area_fractions.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(count_rows[0]));w.writeheader();w.writerows(count_rows)
percent=total_counts/total_counts.sum()*100
with (out/'derived_data'/'phase_area_fractions_pooled.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['class','pixel_count','total_pixels','area_percent'])
 for (name,*_),n,p in zip(classes,total_counts,percent):w.writerow([name,int(n),int(total_counts.sum()),float(p)])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':11,
 'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42})
# Obtain shear compliance from the full matrix, with normal-shear coupling retained.
source=(a.archive/'reference_outputs'/'property_full_homogenized_property_masks_STRICT'/'figure5_property_full_homogenized_property_masks_STRICT_window_values.csv') if a.archive is not None else bundled/'property_window_values.csv'
with source.open() as f: rows=list(csv.DictReader(f))
derived=[]; max_axial_error=0
for r in rows:
 C=np.array([[float(r['C11']),float(r['C12']),float(r['C13'])],
             [float(r['C12']),float(r['C22']),float(r['C23'])],
             [float(r['C13']),float(r['C23']),float(r['C33'])]])
 S=np.linalg.inv(C)
 ex,ey,gxy=1/np.diag(S)
 max_axial_error=max(max_axial_error,abs(ex-float(r['Ex'])),abs(ey-float(r['Ey'])))
 assert gxy>0 and np.isfinite(gxy)
 derived.append({'image':int(r['image']),'window_index':int(r['window_index']),
 'L':int(r['L']),'L_mm':float(r['L_mm']),'Ex':ex,'Ey':ey,'Gxy':gxy})
assert max_axial_error < 1e-8
with (out/'derived_data'/'axial_shear_window_values.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(derived[0]));w.writeheader();w.writerows(derived)
summary=[]
for L in sorted({r['L'] for r in derived}):
 rr=[r for r in derived if r['L']==L]
 row={'L':L,'L_mm':rr[0]['L_mm'],'n_windows':len(rr),'n_images':len({r['image'] for r in rr})}
 for prop in ('Ex','Ey','Gxy'):
  vals=np.array([r[prop] for r in rr]);row[prop+'_mean']=vals.mean();row[prop+'_sem']=vals.std(ddof=1)/np.sqrt(len(vals))
 summary.append(row)
with (out/'derived_data'/'axial_shear_summary.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
# Unified scientific panels: unchanged archived cell values; open isolated
# symbols identify the 384-pixel support represented by image 6 alone.
full_summary=[]
for L in sorted({int(r['L']) for r in rows}):
 rr=[r for r in rows if int(r['L'])==L]
 row={'L':L,'L_mm':float(rr[0]['L_mm']),'n_windows':len(rr),'n_images':len({r['image'] for r in rr})}
 for key in ('K11','K22','C11','C12','C22','C33'):
  vals=np.array([float(r[key]) for r in rr]);row[key+'_mean']=vals.mean();row[key+'_sem']=vals.std(ddof=1)/np.sqrt(len(vals))
 full_summary.append(row)
with (out/'derived_data'/'thermal_stiffness_summary.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(full_summary[0]));w.writeheader();w.writerows(full_summary)

def phase_panel(ax):
 order=np.argsort(-percent); y=np.arange(len(order))
 colors=['#315a7d' if classes[j][0]=='Quartz' else '#8e9ca8' if classes[j][0]=='Pores/background' else '#658f9b' for j in order]
 ax.barh(y,percent[order],height=.66,color=colors)
 ax.set_yticks(y,[classes[j][0] for j in order],fontsize=8)
 ax.invert_yaxis();ax.set_xlim(0,65);ax.set_xticks(np.arange(0,61,10))
 ax.set_xlabel('Area fraction (%)',fontsize=10)
 ax.xaxis.grid(True,alpha=.18);ax.set_axisbelow(True)
 for ypos,j in enumerate(order):ax.text(percent[j]+.6,ypos,f'{percent[j]:.2f}%',va='center',fontsize=8)
 ax.set_title('(a)',loc='left',fontsize=13,fontweight='bold',pad=6)

def property_panel(ax,kind):
 if kind=='b':
  data=full_summary; series=[('K11',r'$\kappa_x^{\rm app}$','o','#1f77b4'),('K22',r'$\kappa_y^{\rm app}$','s','#ff7f0e')]
  ylabel=r'Apparent conductivity (W m$^{-1}$ K$^{-1}$)'
 elif kind=='c':
  data=full_summary; series=[('C11',r'$C_{11}^{\rm app}$','o','#1f77b4'),('C12',r'$C_{12}^{\rm app}$','s','#ff7f0e'),('C22',r'$C_{22}^{\rm app}$','^','#269575'),('C33',r'$C_{33}^{\rm app}$','D','#a465a1')]
  ylabel='Apparent stiffness (GPa)'
 else:
  data=summary; series=[('Ex',r'$E_x^{\rm app}$','o','#1f77b4'),('Ey',r'$E_y^{\rm app}$','s','#ff7f0e'),('Gxy',r'$G_{xy}^{\rm app}$','^','#269575')]
  ylabel='Apparent modulus (GPa)'
 x=np.array([r['L_mm'] for r in data]);same=np.array([r['n_images']==3 for r in data])
 for prop,label,marker,color in series:
  mean=np.array([r[prop+'_mean'] for r in data]);sem=np.array([r[prop+'_sem'] for r in data]);lower=np.minimum(sem,.999*np.maximum(mean,1e-30))
  ax.errorbar(x[same],mean[same],yerr=[lower[same],sem[same]],marker=marker,lw=1.3,ms=3.6,capsize=2,label=label,color=color)
  ax.errorbar(x[~same],mean[~same],yerr=[lower[~same],sem[~same]],marker=marker,linestyle='none',markerfacecolor='white',markeredgewidth=1.1,ms=4.3,capsize=2,color=color)
 ax.axvline(2.04,ls='--',color='black',lw=1.1)
 ax.set_xscale('log');ax.set_yscale('log')
 ax.set_xlabel(r'Window side length, $\ell$ (mm)',fontsize=10);ax.set_ylabel(ylabel,fontsize=10)
 ax.grid(True,which='major',alpha=.24);ax.grid(True,which='minor',alpha=.11)
 ax.legend(fontsize=8.5,loc='upper right',frameon=False)
 ax.set_title(f'({kind})',loc='left',fontsize=13,fontweight='bold',pad=6)
 ax.text(.98,.46,'Open symbols:\nimage 6 only',transform=ax.transAxes,ha='right',va='center',fontsize=7.3,color='#444444')
 ax.tick_params(labelsize=9)

for kind,name in [('a','Figure_5a_linear'),('b','Figure_5b_revised'),('c','Figure_5c_revised'),('d','Figure_5d_axial_shear')]:
 fig,ax=plt.subplots(figsize=(5.4,4.25))
 phase_panel(ax) if kind=='a' else property_panel(ax,kind)
 fig.tight_layout()
 for ext in ('pdf','png'):fig.savefig(out/'IMG'/f'{name}.{ext}',dpi=400,bbox_inches='tight')
 plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(10,7.8))
phase_panel(axes[0,0]); property_panel(axes[0,1],'b');property_panel(axes[1,0],'c');property_panel(axes[1,1],'d')
fig.subplots_adjust(left=.155,right=.99,bottom=.08,top=.96,wspace=.42,hspace=.35)
for ext in ('pdf','png'):fig.savefig(out/'IMG'/f'Figure_3_revised.{ext}',dpi=400,bbox_inches='tight')
plt.close(fig)

# Fig S1: maps reconstructed from already-classified archived numerical arrays.
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
if a.archive is not None:
 k,E,nu=[np.loadtxt(a.archive/'QEMSCAN'/f'{name}Map_1.csv',delimiter=',') for name in ('kappa','E','nu')]
else:
 with np.load(bundled/'qemscan_image1_maps.npz') as maps:k,E,nu=[maps[name] for name in ('kappa','E','nu')]
index=np.full(k.shape,-1,dtype=np.int8)
for j,(_,k0,e0,n0) in enumerate(classes):
 mask=np.isclose(E,e0,atol=1e-8,rtol=0)&np.isclose(nu,n0,atol=1e-8,rtol=0)
 if k0 is not None:mask &= np.isclose(k,k0,atol=1e-8,rtol=0)
 index[mask]=j
assert np.all(index>=0)
rgb=np.array([[255,192,192],[255,0,255],[128,0,0],[255,255,0],[0,255,255],[255,165,0],[0,255,0],[159,94,255],[255,0,0],[255,255,255],[0,0,0]])/255
cmap=ListedColormap(rgb);norm=BoundaryNorm(np.arange(-.5,len(classes)+.5),len(classes))
fig,axes=plt.subplots(1,3,figsize=(11.1,4.15))
extent=(0,k.shape[1]*.01,k.shape[0]*.01,0)
axes[0].imshow(index,cmap=cmap,norm=norm,interpolation='nearest',extent=extent)
axes[0].set_title('(a) Reconstructed class map',fontsize=10,loc='left')
im=axes[1].imshow(k,cmap='viridis',interpolation='nearest',extent=extent,vmin=0,vmax=6.5)
axes[1].set_title('(b) Thermal conductivity',fontsize=10,loc='left')
cb=fig.colorbar(im,ax=axes[1],fraction=.046,pad=.025);cb.set_label(r'W m$^{-1}$ K$^{-1}$',fontsize=8);cb.ax.tick_params(labelsize=8)
im=axes[2].imshow(E,cmap='magma',interpolation='nearest',extent=extent,vmin=0,vmax=100)
axes[2].set_title("(c) Assigned Young's modulus",fontsize=10,loc='left')
cb=fig.colorbar(im,ax=axes[2],fraction=.046,pad=.025);cb.set_label('GPa',fontsize=8);cb.ax.tick_params(labelsize=8)
for ax in axes:
 ax.set_xlabel('x (mm)',fontsize=9);ax.set_ylabel('y (mm)',fontsize=9);ax.tick_params(labelsize=8)
handles=[Patch(facecolor=rgb[j],edgecolor='#777777',linewidth=.5,label=name) for j,(name,*_) in enumerate(classes)]
fig.legend(handles=handles,loc='lower center',ncol=6,fontsize=7.3,frameon=False,bbox_to_anchor=(.5,.0),handlelength=1.1,columnspacing=1)
fig.subplots_adjust(left=.045,right=.97,top=.92,bottom=.23,wspace=.33)
for ext in ('pdf','png'):fig.savefig(out/'IMG'/f'QEMSCAN_property_maps_reconstructed.{ext}',dpi=350,bbox_inches='tight')
plt.close(fig)

print(f'Saved revised figures, {int(total_counts.sum())} classified pixels and {len(rows)} archived cell responses.')
