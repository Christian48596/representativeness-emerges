#!/usr/bin/env python3
"""Rebuild shorter-reference checks and Figure 2 inputs from safe checkpoints.
Added 2026-10-29; original analysis routines are imported without modification.
The 2560-pixel branch explicitly bypasses the production availability gate and
is retained only as an exploratory diagnostic, not a manuscript conclusion.
"""
from pathlib import Path
import argparse
import numpy as np
import rea_validation_wrapper as v

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--checkpoint-dir',type=Path,default=v.CHECKPOINTS)
 p.add_argument('--output-dir',type=Path,default=v.OUT)
 a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
 data=[v.load_checkpoint(a.checkpoint_dir/f'image_{i}.json.gz') for i in range(1,8)]
 extra=[]
 for ref in [1536,1792]:
  sizes=[L for L in v.BASE if L<=ref and all(d['availability49'][L] for d in data)]
  _,_,r=v.summarize_case(data,'masked121',sizes,ref)
  r['enforces_original49_common_size_gate']=True;extra.append(r)
 for ref in v.EXT:
  all_available=all(next(r for r in d['modes']['masked121'] if r['L']==ref)['available'] for d in data)
  if all_available:
   sizes=[L for L in v.BASE+v.EXT if L<=ref and all(next(r for r in d['modes']['masked121'] if r['L']==L)['available'] for d in data)]
   _,_,r=v.summarize_case(data,'masked121',sizes,ref)
   r['enforces_original49_common_size_gate']=False;extra.append(r)
  else:extra.append({'reference':ref,'available_all_images_on121':False})
 v.write_result(a.output_dir/'reference_sensitivity_extra.json',extra)
 per_image,ensemble,base=v.summarize_case(data,'masked121',v.BASE,2048)
 plot={key:np.array([r[key] for r in ensemble]) for key in ['L','sample_mean_Z','eta_C','eta_C_sem']}
 plot.update(tau_Z=v.m.TAU_Z,tau_C=v.m.TAU_C,accepted_L=base['selected'],reference_L=2048)
 for r in v.m.build_ensemble_spectra(per_image,v.BASE):
  if r['L'] in [256,512,1024,2048]:
   n=r['L'];plot[f'k_{n}']=r['spectrum_k'];plot[f'p_{n}']=r['spectrum_p_mean'];plot[f'sem_{n}']=r['spectrum_p_sem']
 np.savez_compressed(a.output_dir/'figure2_plot_data.npz',**plot)
 print('Rebuilt reference_sensitivity_extra.json and figure2_plot_data.npz.')
if __name__=='__main__':main()
