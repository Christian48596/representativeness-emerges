#!/usr/bin/env python3
"""Rebuild reference-sensitivity diagnostics and Figure 2 arrays for N images.

Production results retain the production 49-position availability gate. Larger
reference diagnostics using the sampling grid are explicitly exploratory.
"""
import argparse
import numpy as np
import structural_uncertainty as v


def _available(data, mode, size, gate=True):
    for image in data:
        rows = {r['L']: r for r in image['modes'].get(mode, [])}
        if not rows.get(size, {}).get('available', False):
            return False
        if gate and not image['availability49'].get(size, False):
            return False
    return True


def _case(data, mode, sizes, reference):
    if reference not in sizes or len(sizes) < 4:
        return [], [], dict(
            mode=mode, sizes=sizes, reference=reference, selected=None, table=[],
            reason=('Requested reference is unavailable' if reference not in sizes
                    else 'Fewer than four common valid supports'))
    return v.summarize_case(data, mode, sizes, reference)


def rebuild(data):
    mode = v.baseline_mode()
    extra = []
    configured_short = v.CONFIG['short_reference_sizes']
    short_references = [L for L in configured_short if L in v.BASE and L < v.REFERENCE]
    if not short_references:
        short_references = [L for L in v.BASE if L < v.REFERENCE][-2:]
    for reference in short_references:
        sizes = [L for L in v.BASE if L <= reference and _available(data, mode, L)]
        _, _, result = _case(data, mode, sizes, reference)
        result['enforces_production49_common_size_gate'] = True
        extra.append(result)
    for reference in v.EXT:
        if _available(data, mode, reference, gate=False):
            sizes = sorted({L for L in v.BASE+v.EXT
                            if L <= reference and _available(data, mode, L, gate=False)})
            _, _, result = _case(data, mode, sizes, reference)
            result['enforces_production49_common_size_gate'] = False
            extra.append(result)
        else:
            key = ('available_all_images_on121' if v.CONFIG['n_window_samples'] == 121
                   else 'available_all_images_on_sampling_grid')
            extra.append({'reference': reference, key: False})
    v.write_result(v.OUT/'reference_sensitivity_extra.json', extra)

    sizes = [L for L in v.BASE if _available(data, mode, L)]
    reference = v.REFERENCE if v.REFERENCE in sizes else (sizes[-1] if sizes else v.REFERENCE)
    per_image, ensemble, baseline = _case(data, mode, sizes, reference)
    if ensemble:
        plot = {key: np.array([row[key] for row in ensemble], dtype=float)
                for key in ['L', 'sample_mean_Z', 'eta_C', 'eta_C_sem']}
    else:
        # Keep the requested support axis so a plot can state that every
        # estimate is unavailable. NaNs are markers, not synthetic statistics.
        plot = {'L': np.asarray(v.BASE, dtype=float)}
        for key in ['sample_mean_Z', 'eta_C', 'eta_C_sem']:
            plot[key] = np.full(len(v.BASE), np.nan, dtype=float)
    # Numeric NaN marks no selection. It avoids an object array and preserves
    # allow_pickle=False loading even when no support can be selected.
    plot.update(tau_Z=float(v.m.TAU_Z), tau_C=float(v.m.TAU_C),
                accepted_L=float(baseline['selected']) if baseline['selected'] is not None else np.nan,
                reference_L=int(reference))
    if not v.is_bundled_study(data):
        plot['n_images'] = int(len(data))
        pixel_size = v.CONFIG.get('pixel_size_mm')
        plot['pixel_size_mm'] = float(pixel_size) if pixel_size is not None else np.nan
    if ensemble:
        configured = v.CONFIG['selected_spectra']
        selected = [L for L in configured if L in sizes]
        if not selected:
            indices = np.unique(np.linspace(0, len(sizes)-1, min(4, len(sizes))).round().astype(int))
            selected = [sizes[i] for i in indices]
        for row in v.m.build_ensemble_spectra(per_image, sizes):
            L = row['L']
            if L in selected:
                plot[f'k_{L}'] = row['spectrum_k']
                plot[f'p_{L}'] = row['spectrum_p_mean']
                plot[f'sem_{L}'] = row['spectrum_p_sem']
    np.savez_compressed(v.OUT/'figure2_plot_data.npz', **plot)
    print('Rebuilt reference_sensitivity_extra.json and figure2_plot_data.npz '
          f'for {len(data)} images.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    v.add_common_arguments(parser, actions=False)
    args = parser.parse_args()
    try:
        v.configure(args)
        data = v.load_dataset_checkpoints('image')
        rebuild(data)
    except (ValueError, FileNotFoundError) as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
