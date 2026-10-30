"""Dataset and parameter handling for N-image validation.

No numerical analysis is performed here. A checkpoint manifest prevents an
incomplete or mixed image set being silently treated as a complete ensemble.
"""
from pathlib import Path
import argparse, hashlib, json, re, math
from checkpoint_io import load_checkpoint

DEFAULTS = dict(window_sizes=[256,384,512,640,768,896,1024,1280,1536,1792,2048],
 reference_L=2048,extra_reference_sizes=[2560,3072,3584,4096],
 n_window_samples=121,n_common_k_grid=600,tau_Z=.025,tau_C=.30,
 min_mask_fraction=.98,tail_fraction=.5,low_k_fraction=.25,low_k_weight='inverse',
 screening_sizes=[1536,3072],screening_stride=512,screening_tau_mu=.02,
 screening_tau_sigma=.06,min_component_windows=4,
 short_reference_sizes=[1536,1792],selected_spectra=[256,512,1024,2048],pixel_size_mm=None)

def int_list(text):
    if text is None:return None
    try: values=[int(x.strip()) for x in text.split(',') if x.strip()]
    except ValueError as exc:raise argparse.ArgumentTypeError('Expected comma-separated positive integers.') from exc
    if any(x<=0 for x in values):raise argparse.ArgumentTypeError('Sizes must be positive integers.')
    return sorted(set(values))

def natural_key(path):
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r'(\d+)',str(path))]

def add_options(parser):
    parser.add_argument('--images',nargs='+',type=Path,help='Ordered input images; any positive number and arbitrary filenames.')
    parser.add_argument('--masks',nargs='+',type=Path,help='Ordered masks, one per image.')
    parser.add_argument('--mask-dir',type=Path,help='Mask directory; filenames use image stem plus --mask-suffix.')
    parser.add_argument('--mask-suffix',default='_stationary_mask.tif')
    parser.add_argument('--expected-images',type=int,help='Require this image count when loading checkpoints.')
    for key in ['window_sizes','extra_reference_sizes','screening_sizes','short_reference_sizes','selected_spectra']:
        parser.add_argument('--'+key.replace('_','-'),type=int_list,default=None)
    parser.add_argument('--reference-L',dest='reference_L',type=int,default=None)
    for key in ['n_window_samples','n_common_k_grid','screening_stride','min_component_windows']:
        parser.add_argument('--'+key.replace('_','-'),type=int,default=None)
    for key in ['tau_Z','tau_C','min_mask_fraction','tail_fraction','low_k_fraction','screening_tau_mu','screening_tau_sigma','pixel_size_mm']:
        parser.add_argument('--'+key.replace('_','-'),dest=key,type=float,default=None)
    parser.add_argument('--low-k-weight',choices=['uniform','inverse'],default=None)

def resolved_config(args,base=None):
    result={k:list(v) if isinstance(v,list) else v for k,v in (base or DEFAULTS).items()}
    for key in DEFAULTS:
        value=getattr(args,key,None)
        if value is not None:result[key]=value
    for key,value in result.items():
        if isinstance(value,(int,float)) and not math.isfinite(value):
            raise ValueError(key+' must be finite.')
    if not result['window_sizes'] or len(result['window_sizes'])<4:
        raise ValueError('--window-sizes needs at least four distinct supports.')
    if min(result['window_sizes'])<4:raise ValueError('Support sides must be at least four pixels.')
    if result['reference_L'] not in result['window_sizes']:
        raise ValueError('--reference-L must be included in --window-sizes.')
    for key in ['n_window_samples','screening_stride','min_component_windows']:
        if result[key]<1:raise ValueError(key+' must be positive.')
    if result['n_common_k_grid']<4:raise ValueError('n_common_k_grid must be at least four.')
    for key in ['tau_Z','tau_C']:
        if result[key]<=0:raise ValueError(key+' must be positive.')
    for key in ['screening_tau_mu','screening_tau_sigma']:
        if result[key]<0:raise ValueError(key+' must be nonnegative.')
    for key in ['min_mask_fraction','tail_fraction','low_k_fraction']:
        if not 0<result[key]<=1:raise ValueError(key+' must lie in (0,1].')
    if result['pixel_size_mm'] is not None and result['pixel_size_mm']<=0:
        raise ValueError('--pixel-size-mm must be positive.')
    return result

def image_records(args,root):
    images=getattr(args,'images',None)
    if not images:
        folder=root/'BSE'
        images=sorted((p for p in folder.iterdir() if p.suffix.lower() in {'.tif','.tiff','.png','.jpg','.jpeg','.bmp'}) if folder.is_dir() else [],key=natural_key)
    images=[p.expanduser().resolve() for p in images]
    if not images:raise ValueError('No images found; provide --images PATH... or data-root/BSE/.')
    if len(set(images))!=len(images):raise ValueError('Duplicate image paths are not allowed.')
    masks=getattr(args,'masks',None)
    if masks and getattr(args,'mask_dir',None):raise ValueError('Use --masks or --mask-dir, not both.')
    if masks:
        if len(masks)!=len(images):raise ValueError('--masks must contain one mask for each image.')
        masks=[p.expanduser().resolve() for p in masks]
    else:
        folder=(getattr(args,'mask_dir',None) or root/'reference_outputs/stationary_domains_FULLRES_R2048_STRICT').expanduser().resolve()
        masks=[folder/(p.stem+args.mask_suffix) for p in images]
        if len(set(masks))!=len(masks):raise ValueError('Image stems are duplicated; provide explicit --masks.')
    for p in images:
        if not p.is_file():raise FileNotFoundError(p)
    return [dict(id=i+1,path=p,mask=q) for i,(p,q) in enumerate(zip(images,masks))]

def manifest_for(records,config):
    def identity(path):
        if not path.is_file():return None
        stat=path.stat()
        return {'size_bytes':stat.st_size,'mtime_ns':stat.st_mtime_ns}
    body={'schema_version':1,'images':[{'id':r['id'],'path':str(r['path']),'mask':str(r['mask']),'image_identity':identity(r['path']),'mask_identity':identity(r['mask'])} for r in records],'configuration':config}
    body['fingerprint']=hashlib.sha256(json.dumps(body,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return body

def ensure_manifest(directory,manifest):
    directory.mkdir(parents=True,exist_ok=True)
    path=directory/'dataset_manifest.json'
    if path.exists():
        old=json.loads(path.read_text())
        if old!=manifest:raise ValueError('Checkpoint directory belongs to a different image set or configuration. Use a new --checkpoint-dir.')
    elif any(directory.glob('image_*.json.gz')) or any(directory.glob('rstat_image_*.json.gz')):
        raise ValueError('Checkpoint directory already contains records without a matching dataset manifest. Use a new --checkpoint-dir.')
    else:path.write_text(json.dumps(manifest,indent=2)+'\n')

def discover(directory,prefix,args):
    prefix=prefix.rstrip('_')
    files=sorted(directory.glob(prefix+'_*.json.gz'),key=natural_key)
    if not files:raise ValueError(f'No {prefix} checkpoints found in {directory}.')
    records=[load_checkpoint(p) for p in files]
    ids=[r.get('image') for r in records]
    if any(not isinstance(i,int) or i<=0 for i in ids) or len(set(ids))!=len(ids):
        raise ValueError('Checkpoint image IDs must be distinct positive integers.')
    if ids!=list(range(1,len(ids)+1)):raise ValueError('Checkpoint IDs are incomplete: expected contiguous IDs 1 through N.')
    for path,i in zip(files,ids):
        if path.name!=f'{prefix}_{i}.json.gz':raise ValueError('Checkpoint filename and image ID disagree: '+str(path))
    metadata=[r.get('_dataset_manifest') for r in records]
    if any(m is not None for m in metadata):
        if any(m is None for m in metadata) or any(m!=metadata[0] for m in metadata):
            raise ValueError('Mixed dataset manifests in checkpoints.')
        manifest=metadata[0]
        if len(manifest['images'])!=len(records):raise ValueError('Incomplete checkpoint set: some dataset images have not been processed.')
        configs=[r.get('_configuration') for r in records]
        if any(c!=manifest['configuration'] for c in configs):raise ValueError('Mixed checkpoint configurations.')
        on_disk=directory/'dataset_manifest.json'
        if not on_disk.exists() or json.loads(on_disk.read_text())!=manifest:
            raise ValueError('Missing or mismatched dataset_manifest.json.')
        config=resolved_config(args,manifest['configuration'])
        if config!=manifest['configuration']:raise ValueError('Requested parameters differ from stored checkpoint configuration. Recompute into a new checkpoint directory.')
        if getattr(args,'images',None):
            selected=[str(p.expanduser().resolve()) for p in args.images]
            if selected!=[r['path'] for r in manifest['images']]:raise ValueError('Requested --images do not match checkpoint dataset.')
    else:
        manifest=None
        if getattr(args,'images',None) or getattr(args,'masks',None) or getattr(args,'mask_dir',None):
            raise ValueError('Legacy checkpoints do not identify a custom input list. Recompute images into a new checkpoint directory.')
        config=resolved_config(args)
        if config!=DEFAULTS:raise ValueError('Legacy checkpoints have fixed study parameters. Recompute images before changing configuration.')
    expected=getattr(args,'expected_images',None)
    if expected is not None and expected!=len(records):raise ValueError(f'Expected {expected} images; found {len(records)} checkpoints.')
    return records,config,manifest
