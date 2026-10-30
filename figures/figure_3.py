#!/usr/bin/env python3
"""Assemble manuscript Figure 3: phase composition and apparent properties.

No arguments reproduce the bundled dune-sand study. Explicit map lists or a
--data-root select a new dataset; mineral identities are assigned only when a
--phase-table or --study-preset is explicitly supplied. Numerical cell responses
are read from the solver's per-window CSV; no homogenization is performed here.
"""
from pathlib import Path
import argparse
import csv
import json
import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

REPO = Path(__file__).resolve().parents[1]
BUNDLED = REPO / 'derived_data'
STUDY_CLASSES = [
    ('Quartz', 6.5, 94, .08), ('K-feldspar', 2.3, 70, .25),
    ('Plagioclase', 2., 76, .26), ('Mica', 1., 60, .25),
    ('Calcite', 3.6, 72, .31), ('Dolomite', 5., 94, .28),
    ('Clay minerals', 1.5, 20, .3), ('Gypsum/anhydrite', 2.5, 45, .3),
    ('Heavy minerals', 5., 100, .25), ('Pores/background', .026, .01, .2),
    ('Other/unknown', None, .001, .2)]
STUDY_RGB = np.array([[255,192,192],[255,0,255],[128,0,0],[255,255,0],
    [0,255,255],[255,165,0],[0,255,0],[159,94,255],[255,0,0],
    [255,255,255],[0,0,0]]) / 255
SOURCE_NAME = 'figure5_property_full_homogenized_property_masks_STRICT_window_values.csv'


def natural_key(value):
    return [int(s) if s.isdigit() else s.lower() for s in re.split(r'(\d+)', str(value))]


def load_map(path):
    path = Path(path)
    value = np.load(path, allow_pickle=False) if path.suffix.lower() == '.npy' else np.loadtxt(path, delimiter=',', ndmin=2)
    value = np.asarray(value, dtype=float)
    if value.ndim != 2 or not value.size or not np.all(np.isfinite(value)):
        raise ValueError(f'Map must be a nonempty finite 2D array: {path}')
    return value


def read_maps(paths):
    fields = [load_map(p) for p in paths]
    if len({field.shape for field in fields}) != 1:
        raise ValueError(f'Conductivity, E and nu map shapes differ: {paths}')
    return fields


def discover_maps(root):
    directory = root / 'QEMSCAN'
    groups = []
    for prefix in ['kappaMap_', 'EMap_', 'nuMap_']:
        group = {}
        for pattern in [f'{prefix}*.csv', f'{prefix}*.npy']:
            for path in directory.glob(pattern):
                key = path.stem[len(prefix):]
                if key in group:
                    raise ValueError(f'Ambiguous CSV/NPY map for {prefix}{key}; supply explicit map lists.')
                group[key] = path
        groups.append(group)
    keys = set(groups[0])
    if not keys or any(set(group) != keys for group in groups[1:]):
        raise ValueError('Data root must contain matching QEMSCAN/kappaMap_*, EMap_* and nuMap_* files, or use explicit map lists.')
    ids = sorted(keys, key=natural_key)
    return ids, [[group[key] for group in groups] for key in ids]


def phase_classes(path):
    if path.suffix.lower() == '.json':
        value = json.loads(path.read_text())
        records = value.get('classes', []) if isinstance(value, dict) else value
    else:
        with path.open(newline='') as stream:
            records = list(csv.DictReader(stream))
    result = []
    for row in records:
        if not isinstance(row, dict) or not {'name', 'kappa', 'E', 'nu'}.issubset(row):
            raise ValueError('Each phase-table row needs name, kappa, E and nu fields.')
        vals = []
        for key in ['kappa', 'E', 'nu']:
            val = row[key]
            vals.append(None if val is None or val == '' else float(val))
        if not str(row['name']).strip() or any(v is not None and not np.isfinite(v) for v in vals):
            raise ValueError('Phase names must be nonempty and property matches finite or null.')
        result.append((str(row['name']), *vals))
    if not result or len({r[0] for r in result}) != len(result):
        raise ValueError('Phase table must contain uniquely named classes.')
    return result


def classify(fields, classes, image_id):
    index = np.full(fields[0].shape, -1, dtype=np.int32)
    seen = np.zeros(fields[0].shape, dtype=np.int32)
    for j, (_, *values) in enumerate(classes):
        mask = np.ones(index.shape, dtype=bool)
        for field, value in zip(fields, values):
            if value is not None:
                mask &= np.isclose(field, value, rtol=0, atol=1e-8)
        seen += mask
        index[mask] = j
    if not np.all(seen == 1):
        raise ValueError(f'Phase table leaves {np.count_nonzero(seen == 0)} unmatched and {np.count_nonzero(seen > 1)} ambiguous pixels in image {image_id}.')
    return index


def write_rows(path, rows):
    if not rows:
        return
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sample_sem(values):
    return float(values.std(ddof=1) / np.sqrt(len(values))) if len(values) > 1 else float('nan')


def read_responses(source, pixel_size):
    with source.open(newline='') as stream:
        reader = csv.DictReader(stream)
        required = {'image','window_index','L','K11','K22','C11','C12','C13','C22','C23','C33'}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f'Window-response CSV is missing columns: {sorted(missing)}')
        rows = list(reader)
    if not rows:
        raise ValueError('The window-response CSV contains no solved windows; no property response can be plotted.')
    keys = set()
    lengths = {}
    derived = []
    max_axial_error = 0.
    for r in rows:
        image_id = r['image']
        if not image_id:
            raise ValueError('Each cell response needs a nonempty image identifier.')
        L = int(r['L'])
        if L <= 0:
            raise ValueError('Window sides must be positive.')
        key = (image_id, L, r['window_index'])
        if key in keys:
            raise ValueError(f'Duplicate image/support/window key: {key}')
        keys.add(key)
        physical = r.get('L_mm', '')
        L_mm = float(physical) if physical else (L * pixel_size if pixel_size is not None else None)
        if L_mm is None or not np.isfinite(L_mm) or L_mm <= 0:
            raise ValueError('Supply positive L_mm values or --property-pixel-size-mm for every cell support.')
        if pixel_size is not None and physical and not np.isclose(L_mm, L * pixel_size, rtol=1e-9, atol=1e-12):
            raise ValueError(f'L_mm disagrees with --property-pixel-size-mm for L={L}.')
        if L in lengths and not np.isclose(lengths[L], L_mm, rtol=1e-9, atol=1e-12):
            raise ValueError(f'A common support L={L} has inconsistent physical lengths.')
        lengths[L] = L_mm
        C = np.array([[float(r['C11']),float(r['C12']),float(r['C13'])],
                      [float(r['C12']),float(r['C22']),float(r['C23'])],
                      [float(r['C13']),float(r['C23']),float(r['C33'])]])
        if not np.all(np.isfinite(C)):
            raise ValueError(f'Non-finite stiffness matrix at {key}.')
        try:
            S = np.linalg.inv(C)
        except np.linalg.LinAlgError as exc:
            raise ValueError(f'Singular stiffness matrix at {key}.') from exc
        ex, ey, gxy = 1 / np.diag(S)
        if not np.all(np.isfinite([ex, ey, gxy])) or min(ex, ey, gxy) <= 0:
            raise ValueError(f'Nonpositive or undefined compliance-derived modulus at {key}.')
        for expected, field in [(ex, 'Ex'), (ey, 'Ey')]:
            if r.get(field, ''):
                error = abs(expected - float(r[field]))
                max_axial_error = max(max_axial_error, error)
                if not np.isclose(expected, float(r[field]), rtol=1e-8, atol=1e-8):
                    raise ValueError(f'{field} disagrees with full stiffness inversion at {key}.')
        for field in ['K11', 'K22']:
            if not np.isfinite(float(r[field])) or float(r[field]) <= 0:
                raise ValueError(f'{field} must be finite and positive at {key}.')
        r['L'] = L
        r['L_mm'] = L_mm
        derived.append({'image': image_id, 'window_index': r['window_index'],
                        'L': L, 'L_mm': L_mm, 'Ex': ex, 'Ey': ey, 'Gxy': gxy})
    return rows, derived, max_axial_error


def summarize(rows, quantities):
    output = []
    coverage = {}
    for L in sorted({r['L'] for r in rows}):
        group = [r for r in rows if r['L'] == L]
        image_ids = {r['image'] for r in group}
        coverage[L] = image_ids
        row = {'L': L, 'L_mm': group[0]['L_mm'], 'n_windows': len(group), 'n_images': len(image_ids)}
        for key in quantities:
            values = np.array([float(r[key]) for r in group])
            row[key + '_mean'] = values.mean()
            row[key + '_sem'] = sample_sem(values)
        output.append(row)
    return output, coverage


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', '--archive', dest='archive', type=Path,
        default=Path(os.environ['REA_DATA_ROOT']) if os.environ.get('REA_DATA_ROOT') else None)
    parser.add_argument('--kappa-maps', nargs='+', type=Path)
    parser.add_argument('--E-maps', nargs='+', type=Path)
    parser.add_argument('--nu-maps', nargs='+', type=Path)
    parser.add_argument('--window-values', type=Path, help='Per-window response CSV from the property solver.')
    parser.add_argument('--phase-table', type=Path, help='Optional JSON/CSV classes with name,kappa,E,nu; null/empty property is a wildcard.')
    parser.add_argument('--study-preset', action='store_true', help='Use the published dune-sand phase assignments and nominal 2.04 mm marker.')
    parser.add_argument('--property-pixel-size-mm', type=float, help='Property-map pixel size; required for --rea-pixels. Maps use pixels when omitted.')
    marker = parser.add_mutually_exclusive_group()
    marker.add_argument('--rea-mm', type=float)
    marker.add_argument('--rea-pixels', type=float)
    parser.add_argument('--map-index', type=int, default=1, help='One-based map-list position for the map illustration.')
    parser.add_argument('--output-dir', '--output', dest='output', type=Path, default=REPO/'outputs/figures')
    parser.add_argument('--dpi', type=int, default=400)
    a = parser.parse_args()
    if a.phase_table and a.study_preset:
        parser.error('Use --phase-table or --study-preset, not both.')
    explicit_maps = any([a.kappa_maps, a.E_maps, a.nu_maps])
    if explicit_maps and not all([a.kappa_maps, a.E_maps, a.nu_maps]):
        parser.error('Supply --kappa-maps, --E-maps and --nu-maps together.')
    if a.dpi <= 0 or a.map_index <= 0:
        parser.error('--dpi and --map-index must be positive.')
    if a.property_pixel_size_mm is not None and (not np.isfinite(a.property_pixel_size_mm) or a.property_pixel_size_mm <= 0):
        parser.error('--property-pixel-size-mm must be finite and positive.')
    for name in ['rea_mm', 'rea_pixels']:
        value = getattr(a, name)
        if value is not None and (not np.isfinite(value) or value <= 0):
            parser.error('--rea-mm and --rea-pixels must be finite and positive when supplied.')
    if a.rea_pixels is not None and a.property_pixel_size_mm is None:
        parser.error('--rea-pixels requires --property-pixel-size-mm.')
    # Only the argument-free/compact path automatically uses study assignments.
    compact_study = not explicit_maps and a.archive is None and a.window_values is None and a.phase_table is None
    study = compact_study or a.study_preset
    pixel_size = a.property_pixel_size_mm if a.property_pixel_size_mm is not None else (.01 if study else None)
    rea_mm = a.rea_mm if a.rea_mm is not None else (a.rea_pixels * a.property_pixel_size_mm if a.rea_pixels is not None else (2.04 if study else None))
    classes = phase_classes(a.phase_table) if a.phase_table else (STUDY_CLASSES if study else None)
    map_sets = []
    map_ids = []
    if explicit_maps:
        if len(a.kappa_maps) != len(a.E_maps) or len(a.kappa_maps) != len(a.nu_maps):
            parser.error('The three property-map lists must have the same length and matching order.')
        map_sets = list(zip(a.kappa_maps, a.E_maps, a.nu_maps))
        map_ids = [str(i + 1) for i in range(len(map_sets))]
    elif a.archive is not None:
        map_ids, map_sets = discover_maps(a.archive)
    if map_sets and a.map_index > len(map_sets):
        parser.error(f'--map-index exceeds the {len(map_sets)} supplied map sets.')
    if classes and not map_sets and not compact_study:
        parser.error('Phase classification requires numerical maps; supply map lists or --data-root.')
    source = a.window_values
    if source is None:
        if a.archive is not None:
            source = a.archive/'reference_outputs/property_full_homogenized_property_masks_STRICT'/SOURCE_NAME
        elif compact_study:
            source = BUNDLED/'property_window_values.csv'
        else:
            parser.error('Explicit map inputs require --window-values or --data-root for the solver CSV.')
    out = a.output
    out.mkdir(parents=True, exist_ok=True)
    (out/'derived_data').mkdir(parents=True, exist_ok=True)
    # Remove only this plotter's obsolete class outputs when reusing an output
    # directory for unclassified data; otherwise old mineral labels would remain.
    if classes is None:
        for relative in ['derived_data/phase_area_fractions.csv', 'derived_data/phase_area_fractions_pooled.csv',
                         'Figure_3d.png', 'Figure_3d.pdf']:
            stale = out / relative
            if stale.is_file():
                stale.unlink()
    count_rows = []
    total_counts = np.zeros(len(classes), dtype=np.int64) if classes else None
    display_fields = None
    display_index = None
    if map_sets:
        for i, (image_id, paths) in enumerate(zip(map_ids, map_sets), 1):
            fields = read_maps(paths)
            if classes:
                index = classify(fields, classes, image_id)
                for j, (name, *_) in enumerate(classes):
                    count = int(np.count_nonzero(index == j))
                    total_counts[j] += count
                    count_rows.append({'image':image_id, 'class':name, 'pixel_count':count,
                                       'total_pixels':index.size, 'area_fraction':count/index.size})
            else:
                index = None
            if i == a.map_index:
                display_fields = fields
                display_index = index
    elif compact_study:
        with (BUNDLED/'phase_area_fractions.csv').open() as stream:
            for r in csv.DictReader(stream):
                count_rows.append({'image':int(r['image']), 'class':r['class'],
                    'pixel_count':int(r['pixel_count']), 'total_pixels':int(r['total_pixels']),
                    'area_fraction':float(r['area_fraction'])})
        for j, (name, *_) in enumerate(classes):
            total_counts[j] = sum(r['pixel_count'] for r in count_rows if r['class'] == name)
        with np.load(BUNDLED/'qemscan_image1_maps.npz', allow_pickle=False) as values:
            display_fields = [values[name] for name in ['kappa', 'E', 'nu']]
        display_index = classify(display_fields, classes, 'bundled image 1')
        if a.map_index != 1:
            parser.error('The compact bundle contains maps only for image 1; supply full map inputs for another --map-index.')
    if display_fields is None:
        for extension in ['png', 'pdf']:
            stale = out/f'property_maps.{extension}'
            if stale.is_file():
                stale.unlink()
    percent = total_counts / total_counts.sum() * 100 if classes else None
    if classes:
        write_rows(out/'derived_data/phase_area_fractions.csv', count_rows)
        with (out/'derived_data/phase_area_fractions_pooled.csv').open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['class','pixel_count','total_pixels','area_percent'])
            for (name,*_), count, value in zip(classes, total_counts, percent):
                writer.writerow([name, int(count), int(total_counts.sum()), float(value)])
    rows, derived, max_axial_error = read_responses(source, pixel_size)
    write_rows(out/'derived_data/axial_shear_window_values.csv', derived)
    summary, coverage = summarize(derived, ['Ex','Ey','Gxy'])
    full_summary, _ = summarize(rows, ['K11','K22','C11','C12','C22','C33'])
    write_rows(out/'derived_data/axial_shear_summary.csv', summary)
    write_rows(out/'derived_data/thermal_stiffness_summary.csv', full_summary)
    all_images = set().union(*coverage.values())
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':11,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42})

    def phase_panel(ax):
        order = np.argsort(-percent)
        y = np.arange(len(order))
        colors = ['#315a7d' if classes[j][0]=='Quartz' else '#8e9ca8' if classes[j][0]=='Pores/background' else '#658f9b' for j in order]
        ax.barh(y, percent[order], height=.66, color=colors)
        ax.set_yticks(y, [classes[j][0] for j in order], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlim(0, 65 if study else min(110, max(10, float(percent.max()) * 1.22)))
        if study:
            ax.set_xticks(np.arange(0,61,10))
        ax.set_xlabel('Area fraction (%)', fontsize=10)
        ax.xaxis.grid(True, alpha=.18)
        ax.set_axisbelow(True)
        for ypos, j in zip(y, order):
            ax.text(percent[j]+.6, ypos, f'{percent[j]:.2f}%', va='center', fontsize=8)
        ax.set_title('(a)', loc='left', fontsize=13, fontweight='bold', pad=6)

    def property_panel(ax, kind, label):
        if kind == 'b':
            data = full_summary
            series = [('K11',r'$\kappa_x^{\rm app}$','o','#1f77b4'),('K22',r'$\kappa_y^{\rm app}$','s','#ff7f0e')]
            ylabel = r'Apparent conductivity (W m$^{-1}$ K$^{-1}$)'
        elif kind == 'c':
            data = full_summary
            series = [('C11',r'$C_{11}^{\rm app}$','o','#1f77b4'),('C12',r'$C_{12}^{\rm app}$','s','#ff7f0e'),('C22',r'$C_{22}^{\rm app}$','^','#269575'),('C33',r'$C_{33}^{\rm app}$','D','#a465a1')]
            ylabel = 'Apparent stiffness (GPa)'
        else:
            data = summary
            series = [('Ex',r'$E_x^{\rm app}$','o','#1f77b4'),('Ey',r'$E_y^{\rm app}$','s','#ff7f0e'),('Gxy',r'$G_{xy}^{\rm app}$','^','#269575')]
            ylabel = 'Apparent modulus (GPa)'
        x = np.array([r['L_mm'] for r in data])
        same = np.array([coverage[r['L']] == all_images for r in data])
        negative = any(r[prop+'_mean'] <= 0 for r in data for prop,*_ in series)
        for prop, text, marker_style, color in series:
            mean = np.array([r[prop+'_mean'] for r in data])
            sem = np.array([r[prop+'_sem'] for r in data])
            lower = sem if negative else np.minimum(sem, .999*np.maximum(mean,1e-30))
            # NaN SEM for a single window means no uncertainty estimate, not zero variance.
            for complete, selected in [(True, same), (False, ~same)]:
                if not np.any(selected):
                    continue
                ax.errorbar(x[selected], mean[selected], yerr=[lower[selected],sem[selected]],
                    marker=marker_style, linestyle='-' if complete else 'none', lw=1.3,
                    ms=3.6 if complete else 4.3, capsize=2, color=color,
                    markerfacecolor=color if complete else 'white',
                    label=text if complete or not np.any(same) else None)
        if rea_mm is not None:
            ax.axvline(rea_mm, ls='--', color='black', lw=1.1)
        ax.set_xscale('log')
        if negative:
            magnitudes = [abs(r[prop+'_mean']) for r in data for prop,*_ in series if r[prop+'_mean'] != 0]
            ax.set_yscale('symlog', linthresh=max(min(magnitudes)*.1, 1e-12) if magnitudes else 1e-12)
        else:
            ax.set_yscale('log')
        ax.set_xlabel(r'Window side length, $\ell$ (mm)', fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.grid(True, which='major', alpha=.24)
        ax.grid(True, which='minor', alpha=.11)
        ax.legend(fontsize=8.5, loc='upper right', frameon=False)
        ax.set_title(f'({label})', loc='left', fontsize=13, fontweight='bold', pad=6)
        if np.any(~same):
            text = 'Open symbols:\nreduced image set'
            if study and all_images == {'1','5','6'} and all(coverage[r['L']] == {'6'} for r in data if coverage[r['L']] != all_images):
                text = 'Open symbols:\nimage 6 only'
            ax.text(.98,.46,text,transform=ax.transAxes,ha='right',va='center',fontsize=7.3,color='#444444')
        ax.tick_params(labelsize=9)

    def save(fig, name, dpi=None):
        for ext in ['pdf', 'png']:
            fig.savefig(out/f'{name}.{ext}', dpi=dpi or a.dpi, bbox_inches='tight')
        plt.close(fig)

    panel_specs = [('b','Figure_3b'),('c','Figure_3c'),('d','Figure_3d')]
    if classes:
        fig, ax = plt.subplots(figsize=(5.4,4.25))
        phase_panel(ax)
        fig.tight_layout()
        save(fig,'Figure_3a')
    for kind, name in panel_specs:
        fig, ax = plt.subplots(figsize=(5.4,4.25))
        label = kind if classes else chr(ord('a')+ord(kind)-ord('b'))
        property_panel(ax,kind,label)
        fig.tight_layout()
        save(fig,f'Figure_3{label}')
    if classes:
        fig, axes = plt.subplots(2,2,figsize=(10,7.8))
        phase_panel(axes[0,0])
        property_panel(axes[0,1],'b','b')
        property_panel(axes[1,0],'c','c')
        property_panel(axes[1,1],'d','d')
        fig.subplots_adjust(left=.155,right=.99,bottom=.08,top=.96,wspace=.42,hspace=.35)
    else:
        fig, axes = plt.subplots(1,3,figsize=(14,4.5))
        for ax, (kind,_), label in zip(axes,panel_specs,'abc'):
            property_panel(ax,kind,label)
        fig.tight_layout()
    save(fig,'Figure_3')
    if display_fields is not None:
        k, E, nu = display_fields
        rgb = STUDY_RGB if study else plt.get_cmap('tab20')(np.linspace(0,1,len(classes))) if classes else None
        count = 3
        fig, axes = plt.subplots(1,count,figsize=(11.1,4.15))
        scale = pixel_size if pixel_size is not None else 1.
        unit = 'mm' if pixel_size is not None else 'pixels'
        extent = (0,k.shape[1]*scale,k.shape[0]*scale,0)
        if classes:
            cmap = ListedColormap(rgb)
            norm = BoundaryNorm(np.arange(-.5,len(classes)+.5),len(classes))
            axes[0].imshow(display_index,cmap=cmap,norm=norm,interpolation='nearest',extent=extent)
            axes[0].set_title('(a) Reconstructed class map',fontsize=10,loc='left')
            specs = [(axes[1],k,'viridis','(b) Thermal conductivity',r'W m$^{-1}$ K$^{-1}$'),
                     (axes[2],E,'magma',"(c) Assigned Young's modulus",'GPa')]
        else:
            specs = [(axes[0],k,'viridis','(a) Thermal conductivity',r'W m$^{-1}$ K$^{-1}$'),
                     (axes[1],E,'magma',"(b) Assigned Young's modulus",'GPa'),
                     (axes[2],nu,'cividis','(c) Poisson ratio','Dimensionless')]
        for ax, field, cmap, title, label in specs:
            limits = {'vmin':0, 'vmax':6.5 if field is k else 100} if study and field is not nu else {}
            im = ax.imshow(field,cmap=cmap,interpolation='nearest',extent=extent,**limits)
            ax.set_title(title,fontsize=10,loc='left')
            cb=fig.colorbar(im,ax=ax,fraction=.046,pad=.025)
            cb.set_label(label,fontsize=8)
            cb.ax.tick_params(labelsize=8)
        for ax in axes:
            ax.set_xlabel(f'x ({unit})',fontsize=9)
            ax.set_ylabel(f'y ({unit})',fontsize=9)
            ax.tick_params(labelsize=8)
        if classes:
            handles=[Patch(facecolor=rgb[j],edgecolor='#777777',linewidth=.5,label=name) for j,(name,*_) in enumerate(classes)]
            fig.legend(handles=handles,loc='lower center',ncol=min(6,len(classes)),fontsize=7.3,frameon=False,bbox_to_anchor=(.5,.0),handlelength=1.1,columnspacing=1)
        fig.subplots_adjust(left=.045,right=.97,top=.92,bottom=.23 if classes else .15,wspace=.33)
        save(fig,'property_maps',350 if a.dpi==400 else a.dpi)
    metadata = {
        'n_property_map_sets':len(map_sets) if map_sets else (len({r['image'] for r in count_rows}) if compact_study else 0),
        'n_solved_windows':len(rows), 'response_image_ids':sorted(all_images,key=natural_key),
        'image_ids_by_support':{str(L):sorted(ids,key=natural_key) for L,ids in coverage.items()},
        'phase_classification':bool(classes), 'study_preset':study,
        'phase_table':str(a.phase_table) if a.phase_table else None,
        'rea_marker_mm':rea_mm, 'property_pixel_size_mm':pixel_size,
        'single_window_sem':'undefined (NaN)', 'max_axial_reconstruction_error':max_axial_error,
        'window_values_source':str(source)}
    (out/'derived_data/figure3_metadata.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
    print(f'Saved property figures for {len(all_images)} response images and {len(rows)} solved windows; phase classification {"enabled" if classes else "omitted (no phase table supplied)"}.')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as exc:
        raise SystemExit(f'Error: {exc}') from exc
