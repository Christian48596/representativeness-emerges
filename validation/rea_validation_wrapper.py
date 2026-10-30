#!/usr/bin/env python3
"""Independent audit wrapper. Does not edit or monkeypatch original source.
Core spectral, residual, aggregation and selection routines imported verbatim.
Window enumeration matches source; means streamed to avoid large copied lists.
Image uint8 exact integer levels are converted to float64 before each spectrum;
integer sums for sample means are exact, equivalent to original float64 means.
"""
from pathlib import Path
import importlib.util,sys,time,json,hashlib,argparse,collections,warnings,os
import numpy as np
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
from checkpoint_io import load_checkpoint,save_checkpoint,write_result,finite_json
REPO=Path(__file__).resolve().parents[1]
ROOT=Path(os.environ.get('REA_DATA_ROOT',str(REPO/'data'))).resolve()
OUT=REPO/'outputs/validation'
CHECKPOINTS=REPO/'validation/checkpoints'
SOURCE=REPO/'rea_stationary_mask_cli_v24_general.py'
spec=importlib.util.spec_from_file_location('original_rea',SOURCE)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
BASE=[256,384,512,640,768,896,1024,1280,1536,1792,2048]
EXT=[2560,3072,3584,4096]
MASKDIR=ROOT/'reference_outputs/stationary_domains_FULLRES_R2048_STRICT'
m.N_WINDOW_SAMPLES=121;m.N_COMMON_K_GRID=600;m.MIN_MASK_FRACTION=.98
m.REFERENCE_L=2048;m.WINDOW_SIZES=BASE;m.ORIGINAL_PIXEL_SIZE=1.

def starts(shape,L,n):
 H,W=shape;axis=int(np.ceil(np.sqrt(max(1,n))))
 ys=np.linspace(0,H-L,axis).round().astype(int);xs=np.linspace(0,W-L,axis).round().astype(int)
 ans={(int(y),int(x)) for y in ys for x in xs};ans.add((H//2-L//2,W//2-L//2))
 return sorted(ans)

def one_image(i):
 t0=time.monotonic();path=ROOT/f'BSE/{i}_cut.tif'
 g=np.asarray(Image.open(path).convert('L'),dtype=np.uint8)
 mask=np.asarray(Image.open(MASKDIR/f'{i}_cut_stationary_mask.tif').convert('L'),dtype=np.uint8)>0
 assert mask.shape==g.shape
 print(f'image {i}: {g.shape}, building source integral mask',flush=True)
 ii=m.binary_integral_image(mask)
 material_scale=m.robust_material_scale(g)
 data={'image':i,'shape':g.shape,'mask_fraction':float(mask.mean()),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'modes':{},'availability49':{}}
 for mode,n in [('masked121',121),('masked49',49),('unmasked121',121)]:
  rows=[]
  for L in BASE+EXT:
   ss=starts(g.shape,L,n)
   if mode!='unmasked121':
    ss=[(y,x) for y,x in ss if m.masked_window_fraction(ii,y,x,L)>=.98]
   if mode=='masked49':data['availability49'][L]=bool(ss)
   if not ss:
    rows.append({'L':L,'available':False,'sample_n':0});continue
   if mode=='unmasked121':
    y=g.shape[0]//2-L//2;x=g.shape[1]//2-L//2
   else:
    y,x=min(ss,key=lambda yx:(yx[0]+.5*L-g.shape[0]/2.)**2+(yx[1]+.5*L-g.shape[1]/2.)**2)
   z=np.asarray(g[y:y+L,x:x+L],dtype=float)
   k,p=m.radial_average_spectrum(z,dx=1.)
   eta_b,mean_z,mean_shell=m.boundary_shell_diagnostic(z,alpha=m.ALPHA_SHELL,material_scale=material_scale)
   vals=np.array([float(np.mean(g[yy:yy+L,xx:xx+L])) for yy,xx in ss])
   std=float(np.std(vals,ddof=1)) if len(vals)>1 else 0.
   rows.append(dict(L=L,available=True,mean_Z=mean_z,mean_shell=mean_shell,sample_mean_Z=float(vals.mean()),sample_std_Z=std,sample_sem_Z=std/np.sqrt(len(vals)),sample_n=len(vals),eta_boundary=eta_b,spectrum_k=k,spectrum_p=p,selected_start=(y,x),candidate_starts=ss,sample_values=vals))
   del z
   print(f'image{i} {mode} L{L}: n{len(ss)}, elapsed {time.monotonic()-t0:.1f}s',flush=True)
  data['modes'][mode]=rows
 data['elapsed_seconds']=time.monotonic()-t0
 save_checkpoint(CHECKPOINTS/f'image_{i}.json.gz',data)
 print(f'CHECKPOINT image{i}: {data["elapsed_seconds"]:.1f}s',flush=True)

def assemble(all_data,mode,sizes,ref):
 q=np.linspace(0.,np.pi,m.N_COMMON_K_GRID);weights=m.trapezoid_weights(q);results=[]
 for image in all_data:
  byL={r['L']:r for r in image['modes'][mode]};rref=byL[ref]
  assert rref['available']
  rows=[]
  for L in sizes:
   r=byL[L].copy();assert r['available']
   r['is_reference']=L==ref
   r['eta_C']=np.nan if L==ref else m.low_k_spectral_residual(r['spectrum_k'],r['spectrum_p'],rref['spectrum_k'],rref['spectrum_p'],q,weights)
   r['delta_Z']=r['eta_Z']=abs(r['sample_mean_Z']-rref['sample_mean_Z'])/max(abs(rref['sample_mean_Z']),1e-30)
   r['eta_Z_centered']=abs(r['mean_Z']-rref['mean_Z'])/max(abs(rref['mean_Z']),1e-30)
   rows.append(r)
  results.append(rows)
 return results

def summarize_case(data,mode,sizes,ref):
 p=assemble(data,mode,sizes,ref);a=m.aggregate_ensemble_results(p,sizes,ref)
 selected=m.find_persistent_rea(a,m.TAU_Z,m.TAU_C,m.TAIL_FRACTION)
 keys=['L','mean_Z','sample_mean_Z','sample_std_Z','sample_sem_Z','delta_Z_raw','delta_Z','eta_C','eta_C_std','eta_C_sem','eta_ens','is_reference']
 return p,a,dict(mode=mode,sizes=sizes,reference=ref,selected=selected,table=[{k:r[k] for k in keys} for r in a])

def aggregate(nboot):
 data=[load_checkpoint(CHECKPOINTS/f'image_{i}.json.gz') for i in range(1,8)]
 p,a,base=summarize_case(data,'masked121',BASE,2048)
 archive_path=ROOT/'reference_outputs/combined_REA_final_STRICT_v24/REA_stationary_domains_STRICT_v24_combined_window_statistics.csv'
 if not archive_path.exists():archive_path=REPO/'validation/inputs/original_combined_window_statistics.csv'
 archive=np.genfromtxt(archive_path,delimiter=',',names=True,dtype=None,encoding=None)
 differences={}
 for ak,rk in [('ensemble_sample_mean_Z','sample_mean_Z'),('eta_Z','delta_Z'),('eta_C','eta_C')]:
  differences[ak]=float(np.nanmax(np.abs(archive[ak]-np.array([r[rk] for r in a]))))
 assert differences['ensemble_sample_mean_Z']<1e-9,differences
 assert differences['eta_C']<1e-10,differences
 assert base['selected']==1536,base['selected']
 print('BASELINE EXACT CHECK',differences,flush=True)
 result={'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'baseline_archive_max_absolute_differences':differences,'baseline':base,'cases':[],'mask_fractions':[d['mask_fraction'] for d in data],'elapsed_image_seconds':[d['elapsed_seconds'] for d in data]}
 for mode in ['masked49','unmasked121']:
  _,_,case=summarize_case(data,mode,BASE,2048);result['cases'].append(case)
 for ref in EXT:
  feasible=all(d['availability49'][ref] for d in data)
  if feasible:
   sizes=[L for L in BASE+EXT if L<=ref and all(d['availability49'][L] for d in data)]
   _,_,case=summarize_case(data,'masked121',sizes,ref);result['cases'].append(case)
  else:result['cases'].append({'mode':'masked121','reference':ref,'feasible_under_source_49_candidate_gate':False,'failed_images':[d['image'] for d in data if not d['availability49'][ref]]})
 t=time.monotonic();rng=np.random.default_rng(20261029);draws=rng.integers(0,7,size=(nboot,7));counts=np.array([np.bincount(d,minlength=7) for d in draws]);unique,inverse=np.unique(counts,axis=0,return_inverse=True)
 selected=[]
 for j,c in enumerate(unique):
  ids=np.repeat(np.arange(7),c);aa=m.aggregate_ensemble_results([p[i] for i in ids],BASE,2048)
  s=m.find_persistent_rea(aa,m.TAU_Z,m.TAU_C,m.TAIL_FRACTION);selected.append(np.inf if s is None else float(s))
  if j%200==0:print(f'bootstrap unique {j}/{len(unique)} elapsed{time.monotonic()-t:.1f}s',flush=True)
 selected=np.array(selected)[inverse];finite=selected[np.isfinite(selected)]
 eta=np.array([[r['eta_C'] for r in rows] for rows in p]);bst_eta=np.mean(eta[draws],axis=1)
 bootstrap=dict(n=nboot,seed=20261029,unique_resamples=len(unique),failure_count=int(np.count_nonzero(~np.isfinite(selected))),failure_fraction=float(np.mean(~np.isfinite(selected))),distribution={str(int(k)):int(v) for k,v in collections.Counter(finite).items()},conditional_percentile_interval_pixels=(np.quantile(finite,[.025,.975],method='inverted_cdf').tolist() if len(finite) else [None,None]),unconditional_percentile_interval_pixels=np.quantile(selected,[.025,.975],method='inverted_cdf').tolist(),conditional_mean_pixels=float(finite.mean()) if len(finite) else None,conditional_variance_pixels2=float(finite.var(ddof=1)) if len(finite)>1 else None,conditional_std_pixels=float(finite.std(ddof=1)) if len(finite)>1 else None,pointwise_spectral_percentile_interval=np.nanquantile(bst_eta,[.025,.975],axis=0).T.tolist(),elapsed_seconds=time.monotonic()-t)
 result['bootstrap']=bootstrap
 loo=[]
 for i in range(7):
  aa=m.aggregate_ensemble_results([p[j] for j in range(7) if j!=i],BASE,2048)
  loo.append({'omitted_image':i+1,'selected':m.find_persistent_rea(aa,m.TAU_Z,m.TAU_C,m.TAIL_FRACTION)})
 result['leave_one_image_out']=loo
 for i,rows in enumerate(p):
  np.savetxt(OUT/f'baseline_per_image_{i+1}.csv',[[r['L'],r['sample_mean_Z'],r['sample_n'],r['eta_C'],r['delta_Z']] for r in rows],delimiter=',',header='L,sample_mean_Z,sample_n,eta_C,raw_individual_mean_residual',comments='')
 write_result(OUT/'validation_results.json',result)
 np.savez_compressed(OUT/'bootstrap_draws.npz',draws=draws,selected_pixels=selected,eta_C=bst_eta)
 print(json.dumps(finite_json(result),indent=2,allow_nan=False),flush=True)

def configure(args):
 global ROOT,OUT,CHECKPOINTS,MASKDIR
 ROOT=args.data_root.resolve();OUT=args.output_dir.resolve()
 OUT.mkdir(parents=True,exist_ok=True)
 CHECKPOINTS=(args.checkpoint_dir or (OUT/'checkpoints' if args.image else REPO/'validation/checkpoints')).resolve()
 MASKDIR=ROOT/'reference_outputs/stationary_domains_FULLRES_R2048_STRICT'

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 action=parser.add_mutually_exclusive_group(required=True)
 action.add_argument('--image',type=int,choices=range(1,8),help='Recompute one image from full-resolution data.')
 action.add_argument('--aggregate',action='store_true',help='Recompute bootstrap and sensitivities from safe checkpoints.')
 parser.add_argument('--bootstrap',type=int,default=10000)
 parser.add_argument('--data-root',type=Path,default=ROOT)
 parser.add_argument('--output-dir',type=Path,default=OUT)
 parser.add_argument('--checkpoint-dir',type=Path,default=None,help='Image writes default to output-dir/checkpoints; aggregate reads bundled checkpoints by default.')
 args=parser.parse_args()
 if args.bootstrap < 2:parser.error('--bootstrap must be at least 2')
 configure(args)
 if args.image:one_image(args.image)
 if args.aggregate:aggregate(args.bootstrap)

if __name__=='__main__':main()
