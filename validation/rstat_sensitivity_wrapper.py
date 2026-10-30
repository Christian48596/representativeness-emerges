#!/usr/bin/env python3
"""Screening-window sensitivity for any number of 2D image fields.

The detector and structural numerical routines are imported unchanged. Image
selection and settings use the same interface as rea_validation_wrapper.py.
"""
import argparse
import importlib.util
import json
import time
import numpy as np
from PIL import Image
import rea_validation_wrapper as v

m = v.m
spec = importlib.util.spec_from_file_location(
    'original_detector', v.SOURCE.with_name('detect_stationary_domains_general.py'))
d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d)


def run_image(i):
    record = next((r for r in v.image_records() if r['id'] == i), None)
    if record is None:
        raise ValueError(f'Image index {i} is not in the configured dataset.')
    t = time.monotonic()
    with Image.open(record['path']) as image:
        g = np.asarray(image.convert('L'), dtype=np.uint8)
    material_scale = m.robust_material_scale(g)
    ans = {'image': i, 'shape': g.shape, 'experiments': {}}
    mode = v.baseline_mode()
    stride = v.CONFIG['screening_stride']
    threshold = m.MIN_MASK_FRACTION
    for R in v.CONFIG['screening_sizes']:
        print(f'image{i} R{R} original detector START', flush=True)
        if R > min(g.shape):
            ans['experiments'][R] = {
                'R': R, 'stride': stride, 'mask_fraction': 0.0,
                'coarse_mask_fraction': 0.0, 'mu_ref': None, 'std_ref': None,
                'availability49': {L: False for L in v.BASE},
                'modes': {mode: [{'L': L, 'available': False, 'sample_n': 0}
                                 for L in v.BASE]},
                'reason': 'Screening window exceeds the image dimensions'}
            continue
        cy, cx, means, stds = d.compute_window_stats_fast(g, R, stride)
        raw, mures, stdres, muref, stdref = d.robust_stationary_grid(
            means, stds, tau_mu=v.CONFIG['screening_tau_mu'],
            tau_sigma=v.CONFIG['screening_tau_sigma'])
        clean = d.clean_stationary_grid(
            raw, min_component_windows=v.CONFIG['min_component_windows'],
            do_open_close=True)
        mask = d.grid_to_pixel_mask(clean, g.shape, cy, cx, R)
        print(f'image{i} R{R} maskfraction{mask.mean():.9f} '
              f'elapsed{time.monotonic()-t:.1f}s', flush=True)
        ii = m.binary_integral_image(mask)
        rows, availability49 = [], {}
        for L in v.BASE:
            starts = v.starts(g.shape, L, v.CONFIG['n_window_samples'])
            starts = [(y, x) for y, x in starts
                      if m.masked_window_fraction(ii, y, x, L) >= threshold]
            availability49[L] = any(
                m.masked_window_fraction(ii, y, x, L) >= threshold
                for y, x in v.starts(g.shape, L, 49))
            if not starts:
                rows.append({'L': L, 'available': False, 'sample_n': 0})
                continue
            y, x = min(starts, key=lambda yx:
                       (yx[0]+.5*L-g.shape[0]/2.)**2
                       + (yx[1]+.5*L-g.shape[1]/2.)**2)
            z = np.asarray(g[y:y+L, x:x+L], dtype=float)
            k, p = m.radial_average_spectrum(z, dx=m.ORIGINAL_PIXEL_SIZE)
            etab, meanz, meanshell = m.boundary_shell_diagnostic(
                z, alpha=m.ALPHA_SHELL, material_scale=material_scale)
            values = np.array([float(np.mean(g[yy:yy+L, xx:xx+L]))
                               for yy, xx in starts])
            std = float(np.std(values, ddof=1)) if len(values) > 1 else 0.
            rows.append(dict(
                L=L, available=True, mean_Z=meanz, mean_shell=meanshell,
                sample_mean_Z=float(values.mean()), sample_std_Z=std,
                sample_sem_Z=std/np.sqrt(len(values)), sample_n=len(values),
                eta_boundary=etab, spectrum_k=k, spectrum_p=p,
                selected_start=(y, x)))
            del z
        ans['experiments'][R] = {
            'R': R, 'stride': stride, 'mask_fraction': float(mask.mean()),
            'coarse_mask_fraction': float(clean.mean()),
            'mu_ref': float(muref), 'std_ref': float(stdref),
            'availability49': availability49, 'modes': {mode: rows}}
        Image.fromarray(mask.astype(np.uint8)*255).save(
            v.OUT/f'image{i}_R{R}_mask.tif', compression='tiff_deflate')
        del ii, mask
        print(f'image{i} R{R} COMPLETE elapsed{time.monotonic()-t:.1f}s', flush=True)
    ans['elapsed_seconds'] = time.monotonic()-t
    v.save_checkpoint(v.checkpoint_path('rstat_image', i), v.attach_metadata(ans))


def aggregate():
    data = v.load_dataset_checkpoints('rstat_image')
    out = []
    mode = v.baseline_mode()
    for R in v.CONFIG['screening_sizes']:
        if any(R not in row['experiments'] for row in data):
            raise ValueError(f'Screening window {R} is missing from a checkpoint.')
        inp = [row['experiments'][R] for row in data]
        by_support = [{r['L']: r for r in row['modes'].get(mode, [])} for row in inp]
        eligible = {
            L: [bool(row['availability49'].get(L, False))
                and bool(rows.get(L, {}).get('available', False))
                for row, rows in zip(inp, by_support)]
            for L in v.BASE}
        sizes = [L for L in v.BASE if all(eligible[L])]
        if len(sizes) < 4:
            out.append(dict(
                R=R, stride=v.CONFIG['screening_stride'], available_sizes=sizes,
                selection=None, reason='Fewer than four common valid supports',
                mask_fractions=[row['mask_fraction'] for row in inp],
                failed_images_by_support={
                    L: [row['image'] for row, ok in zip(data, eligible[L]) if not ok]
                    for L in v.BASE}))
            continue
        reference = v.REFERENCE if v.REFERENCE in sizes else sizes[-1]
        _, _, result = v.summarize_case(inp, mode, sizes, reference)
        result.update(R=R, stride=v.CONFIG['screening_stride'],
                      mask_fractions=[row['mask_fraction'] for row in inp],
                      original_requested_reference=v.REFERENCE)
        out.append(result)
    result = {'experiments': out,
              'elapsed_image_seconds': [row['elapsed_seconds'] for row in data]}
    if not v.is_legacy_dataset(data):
        result.update(n_images=len(data), configuration=v.CONFIG)
    v.write_result(v.OUT/'rstat_sensitivity_results.json', result)
    print(json.dumps(v.finite_json(out), indent=2, allow_nan=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    v.add_common_arguments(parser, actions=True)
    args = parser.parse_args()
    try:
        v.configure(args)
        if args.aggregate:
            aggregate()
        elif args.all_images:
            for record in v.image_records():
                run_image(record['id'])
        else:
            run_image(args.image)
    except (ValueError, FileNotFoundError) as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
