#!/usr/bin/env python3
"""Reproduce or apply structural uncertainty analysis to N two-dimensional images.

The spectral, residual, aggregation and selection routines are shared with
structural analysis. Dataset-specific parameters are stored with checkpoints.
"""
from pathlib import Path
import importlib.util,sys,time,json,hashlib,argparse,collections,warnings,os
import numpy as np
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
from checkpoint_io import load_checkpoint,save_checkpoint,write_result,finite_json
import dataset_config as dc
REPO=Path(__file__).resolve().parents[1]
ROOT=Path(os.environ.get('REA_DATA_ROOT',str(REPO/'data'))).resolve()
OUT=REPO/'outputs/validation'
CHECKPOINTS=REPO/'validation/checkpoints'
SOURCE=REPO/'rea_stationary_mask.py'
spec=importlib.util.spec_from_file_location('structural_analysis',SOURCE)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
CONFIG=dc.resolved_config(argparse.Namespace())
BASE=CONFIG['window_sizes'];EXT=CONFIG['extra_reference_sizes'];REFERENCE=CONFIG['reference_L']
ARGS=argparse.Namespace();RECORDS=[];MANIFEST=None
MASKDIR=ROOT/'reference_outputs/stationary_domains_FULLRES_R2048_STRICT'

def apply_config(config):
 global CONFIG,BASE,EXT,REFERENCE
 CONFIG=config;BASE=config['window_sizes'];EXT=config['extra_reference_sizes'];REFERENCE=config['reference_L']
 m.N_WINDOW_SAMPLES=config['n_window_samples'];m.N_COMMON_K_GRID=config['n_common_k_grid']
 m.MIN_MASK_FRACTION=config['min_mask_fraction'];m.REFERENCE_L=REFERENCE;m.WINDOW_SIZES=BASE
 m.ORIGINAL_PIXEL_SIZE=1.;m.TAU_Z=config['tau_Z'];m.TAU_C=config['tau_C']
 m.TAIL_FRACTION=config['tail_fraction'];m.LOW_K_FRACTION=config['low_k_fraction'];m.LOW_K_WEIGHT=config['low_k_weight']
apply_config(CONFIG)

def baseline_mode():return f"masked{CONFIG['n_window_samples']}"
def image_records():return RECORDS

def checkpoint_path(prefix,i):return CHECKPOINTS/f'{prefix.rstrip("_")}_{i}.json.gz'

def attach_metadata(data):
 data['_configuration']=CONFIG
 data['_dataset_manifest']=MANIFEST
 return data

def is_bundled_study(data):return all('_configuration' not in row for row in data)

def load_dataset_checkpoints(prefix='image'):
 global MANIFEST
 data,config,manifest=dc.discover(CHECKPOINTS,prefix,ARGS)
 apply_config(config);MANIFEST=manifest
 return data

def add_common_arguments(parser,actions=True):
 if actions:
  action=parser.add_mutually_exclusive_group(required=True)
  action.add_argument('--image',type=int,help='Process a 1-based position in the supplied/discovered image list.')
  action.add_argument('--all-images',action='store_true',help='Process every image sequentially.')
  action.add_argument('--aggregate',action='store_true',help='Aggregate the complete checkpoint set.')
 parser.add_argument('--data-root',type=Path,default=ROOT)
 parser.add_argument('--output-dir',type=Path,default=OUT)
 parser.add_argument('--checkpoint-dir',type=Path,default=None)
 dc.add_options(parser)

def configure(args):
 global ROOT,OUT,CHECKPOINTS,MASKDIR,ARGS,RECORDS,MANIFEST
 ARGS=args;ROOT=args.data_root.expanduser().resolve();OUT=args.output_dir.expanduser().resolve()
 OUT.mkdir(parents=True,exist_ok=True)
 writing=bool(getattr(args,'image',None) is not None or getattr(args,'all_images',False))
 CHECKPOINTS=(args.checkpoint_dir or (OUT/'checkpoints' if writing else REPO/'validation/checkpoints')).expanduser().resolve()
 MASKDIR=(getattr(args,'mask_dir',None) or ROOT/'reference_outputs/stationary_domains_FULLRES_R2048_STRICT').expanduser().resolve()
 apply_config(dc.resolved_config(args))
 RECORDS=[];MANIFEST=None
 if writing:
  RECORDS=dc.image_records(args,ROOT)
  if getattr(args,'expected_images',None) is not None and args.expected_images!=len(RECORDS):
   raise ValueError(f'Expected {args.expected_images} input images; found {len(RECORDS)}.')
  i=getattr(args,'image',None)
  if i is not None and not 1<=i<=len(RECORDS):raise ValueError(f'--image must be between 1 and {len(RECORDS)}.')
  MANIFEST=dc.manifest_for(RECORDS,CONFIG);dc.ensure_manifest(CHECKPOINTS,MANIFEST)

def starts(shape,L,n):
 H,W=shape
 if L>min(H,W) or L<1:return []
 axis=int(np.ceil(np.sqrt(max(1,n))))
 ys=np.linspace(0,H-L,axis).round().astype(int);xs=np.linspace(0,W-L,axis).round().astype(int)
 ans={(int(y),int(x)) for y in ys for x in xs};ans.add((H//2-L//2,W//2-L//2))
 return sorted(ans)

def one_image(i):
 record=RECORDS[i-1];t0=time.monotonic();path=record['path']
 g=np.asarray(Image.open(path).convert('L'),dtype=np.uint8)
 if not record['mask'].is_file():raise FileNotFoundError(f"Mask for image {i} not found: {record['mask']}")
 mask=np.asarray(Image.open(record['mask']).convert('L'),dtype=np.uint8)>0
 if mask.shape!=g.shape:raise ValueError(f'Mask shape differs from image {i}.')
 print(f'image {i}: {path.name}, {g.shape}',flush=True)
 ii=m.binary_integral_image(mask);material_scale=m.robust_material_scale(g)
 data={'image':i,'shape':g.shape,'mask_fraction':float(mask.mean()),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'modes':{},'availability49':{}}
 sizes=sorted(set(BASE+EXT));n_samples=CONFIG['n_window_samples']
 for L in sizes:
  data['availability49'][L]=any(m.masked_window_fraction(ii,y,x,L)>=m.MIN_MASK_FRACTION for y,x in starts(g.shape,L,49))
 modes=list(dict.fromkeys([(baseline_mode(),n_samples),('masked49',49),(f'unmasked{n_samples}',n_samples)]))
 for mode,n in modes:
  rows=[]
  for L in sizes:
   ss=starts(g.shape,L,n)
   if not mode.startswith('unmasked'):
    ss=[(y,x) for y,x in ss if m.masked_window_fraction(ii,y,x,L)>=m.MIN_MASK_FRACTION]
   if not ss:rows.append({'L':L,'available':False,'sample_n':0});continue
   if mode.startswith('unmasked'):y=g.shape[0]//2-L//2;x=g.shape[1]//2-L//2
   else:y,x=min(ss,key=lambda yx:(yx[0]+.5*L-g.shape[0]/2.)**2+(yx[1]+.5*L-g.shape[1]/2.)**2)
   z=np.asarray(g[y:y+L,x:x+L],dtype=float);k,p=m.radial_average_spectrum(z,dx=1.)
   eta_b,mean_z,mean_shell=m.boundary_shell_diagnostic(z,alpha=m.ALPHA_SHELL,material_scale=material_scale)
   vals=np.array([float(np.mean(g[yy:yy+L,xx:xx+L])) for yy,xx in ss]);std=float(np.std(vals,ddof=1)) if len(vals)>1 else 0.
   rows.append(dict(L=L,available=True,mean_Z=mean_z,mean_shell=mean_shell,sample_mean_Z=float(vals.mean()),sample_std_Z=std,sample_sem_Z=std/np.sqrt(len(vals)),sample_n=len(vals),eta_boundary=eta_b,spectrum_k=k,spectrum_p=p,selected_start=(y,x),candidate_starts=ss,sample_values=vals))
   del z
  data['modes'][mode]=rows
 data['elapsed_seconds']=time.monotonic()-t0
 save_checkpoint(checkpoint_path('image',i),attach_metadata(data))
 print(f'CHECKPOINT image{i}: {data["elapsed_seconds"]:.1f}s',flush=True)

def assemble(all_data,mode,sizes,ref):
 q=np.linspace(0.,np.pi,m.N_COMMON_K_GRID);weights=m.trapezoid_weights(q);results=[]
 for image in all_data:
  byL={r['L']:r for r in image['modes'][mode]};rref=byL[ref]
  rows=[]
  for L in sizes:
   r=byL[L].copy()
   r['is_reference']=L==ref
   r['eta_C']=np.nan if L==ref else m.low_k_spectral_residual(r['spectrum_k'],r['spectrum_p'],rref['spectrum_k'],rref['spectrum_p'],q,weights)
   r['delta_Z']=r['eta_Z']=abs(r['sample_mean_Z']-rref['sample_mean_Z'])/max(abs(rref['sample_mean_Z']),1e-30)
   r['eta_Z_centered']=abs(r['mean_Z']-rref['mean_Z'])/max(abs(rref['mean_Z']),1e-30)
   rows.append(r)
  results.append(rows)
 return results

def available(data,mode,L):
 return all(any(r['L']==L and r.get('available',False) for r in image['modes'].get(mode,[])) for image in data)

def eligible_sizes(data,mode,sizes,enforce_gate=True):
 return [L for L in sorted(set(sizes)) if available(data,mode,L) and (not enforce_gate or all(row.get('availability49',{}).get(L,False) for row in data))]

def summarize_case(data,mode,sizes,ref):
 summary=dict(mode=mode,sizes=sizes,reference=ref,selected=None,table=[])
 if len(sizes)<4 or ref not in sizes or not data or not all(available(data,mode,L) for L in sizes):
  summary['reason']='Fewer than four common valid supports or reference unavailable.'
  return [],[],summary
 p=assemble(data,mode,sizes,ref);a=m.aggregate_ensemble_results(p,sizes,ref)
 summary['selected']=m.find_persistent_rea(a,m.TAU_Z,m.TAU_C,m.TAIL_FRACTION)
 keys=['L','mean_Z','sample_mean_Z','sample_std_Z','sample_sem_Z','delta_Z_raw','delta_Z','eta_C','eta_C_std','eta_C_sem','eta_ens','is_reference']
 summary['table']=[{k:r[k] for k in keys} for r in a]
 return p,a,summary

def baseline_case(data):
 sizes=eligible_sizes(data,baseline_mode(),BASE)
 ref=REFERENCE if REFERENCE in sizes else (sizes[-1] if sizes else REFERENCE)
 p,a,summary=summarize_case(data,baseline_mode(),sizes,ref)
 if ref!=REFERENCE:summary['requested_reference']=REFERENCE
 return p,a,summary

def verify_reference(data,base,a):
 if len(data)!=7 or CONFIG!=dc.DEFAULTS:raise ValueError('--verify-reference applies only to the recorded seven-image study configuration.')
 # The study check is independent of a user's current external data root.
 path=REPO/'validation/inputs/study_window_statistics.csv'
 reference_table=np.genfromtxt(path,delimiter=',',names=True,dtype=None,encoding=None)
 differences={}
 for ak,rk in [('ensemble_sample_mean_Z','sample_mean_Z'),('eta_Z','delta_Z'),('eta_C','eta_C')]:
  if len(reference_table)!=len(a):raise ValueError('Reference table length differs from calculated support sequence.')
  differences[ak]=float(np.nanmax(np.abs(reference_table[ak]-np.array([r[rk] for r in a]))))
 if differences['ensemble_sample_mean_Z']>=1e-9 or differences['eta_C']>=1e-10 or base['selected']!=1536:
  raise ValueError('Study reference verification failed: '+str(differences))
 return differences

def aggregate(nboot,seed=20261029,verify=False):
 data=load_dataset_checkpoints();bundled_study=is_bundled_study(data);n=len(data)
 p,a,base=baseline_case(data);sizes=base['sizes'];ref=base['reference']
 differences=verify_reference(data,base,a) if verify else {}
 result={'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'baseline_reference_max_absolute_differences':differences,'baseline':base,'cases':[],'mask_fractions':[d['mask_fraction'] for d in data],'elapsed_image_seconds':[d['elapsed_seconds'] for d in data]}
 if not bundled_study:
  result.update(n_images=n,configuration=CONFIG,image_paths=[r['path'] for r in MANIFEST['images']],dataset_fingerprint=MANIFEST['fingerprint'])
 for mode in ['masked49',f"unmasked{CONFIG['n_window_samples']}"]:
  ss=eligible_sizes(data,mode,BASE,enforce_gate=not mode.startswith('unmasked'))
  rr=REFERENCE if REFERENCE in ss else (ss[-1] if ss else REFERENCE)
  _,_,case=summarize_case(data,mode,ss,rr);result['cases'].append(case)
 for reference in EXT:
  feasible=all(d['availability49'].get(reference,False) for d in data)
  if feasible:
   ss=eligible_sizes(data,baseline_mode(),[L for L in sorted(set(BASE+EXT)) if L<=reference])
   _,_,case=summarize_case(data,baseline_mode(),ss,reference);result['cases'].append(case)
  else:result['cases'].append({'mode':baseline_mode(),'reference':reference,'feasible_under_source_49_candidate_gate':False,'failed_images':[d['image'] for d in data if not d['availability49'].get(reference,False)]})
 t=time.monotonic();rng=np.random.default_rng(seed);draws=rng.integers(0,n,size=(nboot,n));counts=np.array([np.bincount(d,minlength=n) for d in draws]);unique,inverse=np.unique(counts,axis=0,return_inverse=True)
 selected=[]
 for j,c in enumerate(unique):
  if p:
   ids=np.repeat(np.arange(n),c);aa=m.aggregate_ensemble_results([p[i] for i in ids],sizes,ref)
   s=m.find_persistent_rea(aa,m.TAU_Z,m.TAU_C,m.TAIL_FRACTION)
  else:s=None
  selected.append(np.inf if s is None else float(s))
  if j%200==0:print(f'bootstrap unique {j}/{len(unique)}',flush=True)
 selected=np.array(selected)[inverse];finite=selected[np.isfinite(selected)]
 eta=np.array([[r['eta_C'] for r in rows] for rows in p]) if p else np.empty((n,0))
 bst_eta=np.mean(eta[draws],axis=1)
 with warnings.catch_warnings():
  warnings.simplefilter('ignore',RuntimeWarning)
  pointwise=np.nanquantile(bst_eta,[.025,.975],axis=0).T.tolist() if bst_eta.shape[1] else []
 bootstrap=dict(n=nboot,seed=seed,unique_resamples=len(unique),failure_count=int(np.count_nonzero(~np.isfinite(selected))),failure_fraction=float(np.mean(~np.isfinite(selected))),distribution={str(int(k)):int(v) for k,v in collections.Counter(finite).items()},conditional_percentile_interval_pixels=(np.quantile(finite,[.025,.975],method='inverted_cdf').tolist() if len(finite) else [None,None]),unconditional_percentile_interval_pixels=np.quantile(selected,[.025,.975],method='inverted_cdf').tolist(),conditional_mean_pixels=float(finite.mean()) if len(finite) else None,conditional_variance_pixels2=float(finite.var(ddof=1)) if len(finite)>1 else None,conditional_std_pixels=float(finite.std(ddof=1)) if len(finite)>1 else None,pointwise_spectral_percentile_interval=pointwise,elapsed_seconds=time.monotonic()-t)
 if not bundled_study:bootstrap.update(n_images=n,independent_image_uncertainty_estimable=n>1)
 if n==1:bootstrap['note']='One image: resampling is degenerate and cannot estimate independent-image uncertainty.'
 result['bootstrap']=bootstrap
 loo=[]
 if n>1:
  for i in range(n):
   if p:
    aa=m.aggregate_ensemble_results([p[j] for j in range(n) if j!=i],sizes,ref);s=m.find_persistent_rea(aa,m.TAU_Z,m.TAU_C,m.TAIL_FRACTION)
   else:s=None
   loo.append({'omitted_image':data[i]['image'],'selected':s})
 result['leave_one_image_out']=loo
 for image,rows in zip(data,p):
  np.savetxt(OUT/f'baseline_per_image_{image["image"]}.csv',[[r['L'],r['sample_mean_Z'],r['sample_n'],r['eta_C'],r['delta_Z']] for r in rows],delimiter=',',header='L,sample_mean_Z,sample_n,eta_C,raw_individual_mean_residual',comments='')
 write_result(OUT/'validation_results.json',result)
 np.savez_compressed(OUT/'bootstrap_draws.npz',draws=draws,selected_pixels=selected,eta_C=bst_eta)
 print(f'Images: {n}; baseline: {base["selected"]}; bootstrap non-selections: {bootstrap["failure_count"]}/{nboot}',flush=True)

def main():
 parser=argparse.ArgumentParser(description=__doc__);add_common_arguments(parser)
 parser.add_argument('--bootstrap',type=int,default=10000)
 parser.add_argument('--seed',type=int,default=20261029)
 parser.add_argument('--verify-reference',action='store_true',help='Verify the recorded seven-image study; not for new datasets.')
 args=parser.parse_args()
 if args.bootstrap<2:parser.error('--bootstrap must be at least 2')
 try:
  configure(args)
  if args.all_images:
   for record in image_records():one_image(record['id'])
  elif args.image is not None:one_image(args.image)
  else:aggregate(args.bootstrap,args.seed,args.verify_reference)
 except (ValueError,FileNotFoundError) as exc:parser.exit(2,f'Error: {exc}\n')
if __name__=='__main__':main()
