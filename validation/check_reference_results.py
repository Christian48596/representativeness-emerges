#!/usr/bin/env python3
"""Recompute the compact revision results and compare all scientific values.
No raw images or FEM solves are needed. About 15 seconds on the audit machine.
"""
from pathlib import Path
import argparse,json,subprocess,sys,tempfile,hashlib
import numpy as np

REPO=Path(__file__).resolve().parents[1]
EXPECTED_GIT_BLOBS={
 'detect_stationary_domains_general.py':'9682ff8bec2d55e2f2481fb361326ed47ee8c943',
 'rea_stationary_mask_cli_v24_general.py':'862983356fc38aa8df56a1a3dea1fc962c150537',
 'property_full_homogenized_inside_property_masks.py':'ca4ae44afa2e691255608a4d62e38dbd3ef92997'}
# Historical source_sha256 hashes include the acquisition copy's extra trailing
# newline. Git originals are restored exactly; scientific values are unaffected.
SKIP={'elapsed_seconds','elapsed_image_seconds','source_sha256'}

def compare(a,b,path='root'):
 if isinstance(a,dict):
  assert set(a)==set(b),(path,set(a)^set(b))
  for k in a:
   if k not in SKIP:compare(a[k],b[k],path+'/'+k)
 elif isinstance(a,list):
  assert len(a)==len(b),path
  for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'/{i}')
 elif isinstance(a,(int,float)) and not isinstance(a,bool):
  assert np.isclose(a,b,rtol=1e-12,atol=1e-12,equal_nan=True),(path,a,b)
 else:assert a==b,(path,a,b)

def main():
 p=argparse.ArgumentParser(description=__doc__);p.parse_args()
 for name,expected in EXPECTED_GIT_BLOBS.items():
  data=(REPO/name).read_bytes();actual=hashlib.sha1(f'blob {len(data)}\0'.encode()+data).hexdigest()
  assert actual==expected,(name,actual)
 with tempfile.TemporaryDirectory(prefix='rea-verification-') as td:
  out=Path(td)
  for script,args in [('rea_validation_wrapper.py',['--aggregate','--bootstrap','10000']),('rstat_sensitivity_wrapper.py',['--aggregate']),('rebuild_derived.py',[])]:
   result=subprocess.run([sys.executable,str(REPO/'validation'/script),*args,'--output-dir',str(out)],capture_output=True,text=True)
   if result.returncode:
    print(result.stdout);print(result.stderr,file=sys.stderr);raise SystemExit(result.returncode)
  for name in ['validation_results.json','rstat_sensitivity_results.json','reference_sensitivity_extra.json']:
   compare(json.loads((REPO/'validation/results'/name).read_text()),json.loads((out/name).read_text()))
  for actual,expected in [(out/'bootstrap_draws.npz',REPO/'validation/results/bootstrap_draws.npz'),(out/'figure2_plot_data.npz',REPO/'derived_data/figure2_plot_data.npz')]:
   with np.load(actual,allow_pickle=False) as a,np.load(expected,allow_pickle=False) as b:
    assert set(a.files)==set(b.files)
    for key in a.files:np.testing.assert_allclose(a[key],b[key],rtol=1e-12,atol=1e-12,equal_nan=True,err_msg=key)
  for i in range(1,8):
   name=f'baseline_per_image_{i}.csv'
   np.testing.assert_allclose(np.loadtxt(out/name,delimiter=',',skiprows=1),np.loadtxt(REPO/'validation/results'/name,delimiter=',',skiprows=1),rtol=1e-12,atol=1e-12,equal_nan=True)
 print('PASS: three original Git blobs; baseline; all 10,000 bootstrap draws; leave-one-out; grid, filtering, reference and R_stat sensitivities; seven per-image CSVs; Figure 2 arrays.')
 print('Baseline: 1536 pixels. Bootstrap non-selections: 3310/10000. Conditional interval: 896–1536 pixels.')
if __name__=='__main__':main()
