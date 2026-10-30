#!/usr/bin/env python3
"""Parameter sensitivity using unchanged original detector and REA functions."""
from pathlib import Path
import importlib.util,time,json,sys,argparse
import numpy as np
from PIL import Image
import rea_validation_wrapper as v
m=v.m
spec=importlib.util.spec_from_file_location('original_detector',v.SOURCE.with_name('detect_stationary_domains_general.py'))
d=importlib.util.module_from_spec(spec);spec.loader.exec_module(d)

def run_image(i):
 t=time.monotonic();g=np.asarray(Image.open(v.ROOT/f'BSE/{i}_cut.tif').convert('L'),dtype=np.uint8)
 material_scale=m.robust_material_scale(g);ans={'image':i,'shape':g.shape,'experiments':{}}
 for R in [1536,3072]:
  print(f'image{i} R{R} original detector START',flush=True)
  cy,cx,means,stds=d.compute_window_stats_fast(g,R,512)
  raw,mures,stdres,muref,stdref=d.robust_stationary_grid(means,stds,tau_mu=.02,tau_sigma=.06)
  clean=d.clean_stationary_grid(raw,min_component_windows=4,do_open_close=True)
  mask=d.grid_to_pixel_mask(clean,g.shape,cy,cx,R)
  print(f'image{i} R{R} maskfraction{mask.mean():.9f} elapsed{time.monotonic()-t:.1f}s',flush=True)
  ii=m.binary_integral_image(mask);rows=[];availability49={}
  for L in v.BASE:
   ss=v.starts(g.shape,L,121);ss=[(y,x) for y,x in ss if m.masked_window_fraction(ii,y,x,L)>=.98]
   availability49[L]=any(m.masked_window_fraction(ii,y,x,L)>=.98 for y,x in v.starts(g.shape,L,49))
   if not ss:rows.append({'L':L,'available':False,'sample_n':0});continue
   y,x=min(ss,key=lambda yx:(yx[0]+.5*L-g.shape[0]/2.)**2+(yx[1]+.5*L-g.shape[1]/2.)**2)
   z=np.asarray(g[y:y+L,x:x+L],dtype=float);k,p=m.radial_average_spectrum(z,dx=1.)
   etab,meanz,meanshell=m.boundary_shell_diagnostic(z,alpha=m.ALPHA_SHELL,material_scale=material_scale)
   vals=np.array([float(np.mean(g[yy:yy+L,xx:xx+L])) for yy,xx in ss]);std=float(np.std(vals,ddof=1)) if len(vals)>1 else 0.
   rows.append(dict(L=L,available=True,mean_Z=meanz,mean_shell=meanshell,sample_mean_Z=float(vals.mean()),sample_std_Z=std,sample_sem_Z=std/np.sqrt(len(vals)),sample_n=len(vals),eta_boundary=etab,spectrum_k=k,spectrum_p=p,selected_start=(y,x)))
   del z
  ans['experiments'][R]={'R':R,'stride':512,'mask_fraction':float(mask.mean()),'coarse_mask_fraction':float(clean.mean()),'mu_ref':float(muref),'std_ref':float(stdref),'availability49':availability49,'modes':{'masked121':rows}}
  # Derived masks compressed to permit review/reproduction without touching originals.
  Image.fromarray(mask.astype(np.uint8)*255).save(v.OUT/f'image{i}_R{R}_mask.tif',compression='tiff_deflate')
  del ii,mask
  print(f'image{i} R{R} COMPLETE elapsed{time.monotonic()-t:.1f}s',flush=True)
 ans['elapsed_seconds']=time.monotonic()-t
 v.save_checkpoint(v.CHECKPOINTS/f'rstat_image_{i}.json.gz',ans)

def aggregate():
 data=[v.load_checkpoint(v.CHECKPOINTS/f'rstat_image_{i}.json.gz') for i in range(1,8)];out=[]
 for R in [1536,3072]:
  inp=[row['experiments'][R] for row in data]
  sizes=[L for L in v.BASE if all(row['availability49'][L] for row in inp)]
  if len(sizes)<4:
   out.append(dict(R=R,stride=512,available_sizes=sizes,selection=None,reason='Fewer than four common valid supports',mask_fractions=[row['mask_fraction'] for row in inp],failed_images_by_support={L:[j+1 for j,row in enumerate(inp) if not row['availability49'][L]] for L in v.BASE}));continue
  ref=2048 if 2048 in sizes else sizes[-1]
  _,_,result=v.summarize_case(inp,'masked121',sizes,ref)
  result.update(R=R,stride=512,mask_fractions=[row['mask_fraction'] for row in inp],original_requested_reference=2048)
  out.append(result)
 v.write_result(v.OUT/'rstat_sensitivity_results.json',{'experiments':out,'elapsed_image_seconds':[r['elapsed_seconds'] for r in data]})
 print(json.dumps(v.finite_json(out),indent=2,allow_nan=False),flush=True)
def main():
 parser=argparse.ArgumentParser(description=__doc__)
 action=parser.add_mutually_exclusive_group(required=True)
 action.add_argument('--image',type=int,choices=range(1,8))
 action.add_argument('--aggregate',action='store_true')
 parser.add_argument('--data-root',type=Path,default=v.ROOT)
 parser.add_argument('--output-dir',type=Path,default=v.OUT)
 parser.add_argument('--checkpoint-dir',type=Path,default=None)
 args=parser.parse_args();v.configure(args)
 if args.image:run_image(args.image)
 else:aggregate()

if __name__=='__main__':main()
