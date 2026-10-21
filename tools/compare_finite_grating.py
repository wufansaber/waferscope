"""Compare finite-grating runs without fitting their amplitude or phase."""
import argparse
import json
from pathlib import Path

import numpy as np

from fullwave_prototype.grating_analysis import relative


def compare(base, refined):
    reports = [json.loads((folder/'report.json').read_text(encoding='utf-8')) for folder in (base, refined)]
    a, b = reports
    if a['geometry'] != b['geometry']:
        raise ValueError('Identical geometry required for a convergence comparison')
    keys = ('wavelength_nm', 'substrate_n', 'substrate_k', 'incidence_deg', 'polarization', 'beam_profile')
    if any(a['parameters'][key] != b['parameters'][key] for key in keys):
        raise ValueError('Identical material and illumination required')
    arrays = []
    for folder in (base, refined):
        with np.load(folder/'cap_checks.npz', allow_pickle=False) as data:
            orders = [int(name[len('volume_q'):]) for name in data.files if name.startswith('volume_q')]
            arrays.append((data['cap_directions'], data['volume_q{}'.format(max(orders))]))
    if not np.allclose(arrays[0][0], arrays[1][0]):
        raise ValueError('Identical observation rays required')
    same_mesh_controls = all(a['parameters'][key] == b['parameters'][key] for key in (
        'air_mesh_um', 'substrate_mesh_um', 'pml_thickness_um', 'pml_mesh_scale'))
    result = dict(geometry=a['geometry'], same_mesh_controls=same_mesh_controls,
        element_counts=[r['solver']['elements'] for r in reports],
        complex_amplitude_relative_l2=relative(arrays[0][1], arrays[1][1]),
        selected_cap_power_change_relative=float(a['checks']['CH3_NA_0p6_DSC_um2']/
            b['checks']['CH3_NA_0p6_DSC_um2']-1),
        status='NOT_CERTIFIED',
        caveat='Differences are diagnostics, not a certified error bound. No amplitude/phase fit; no forced symmetry.',
        runs=[])
    for folder, report in zip((base, refined), reports):
        resources_path = folder.with_suffix('.resources.json')
        result['runs'].append(dict(path=str(folder), parameters=report['parameters'],
            solver=report['solver'], checks=report['checks'],
            resources=json.loads(resources_path.read_text(encoding='utf-8')) if resources_path.exists() else None))
    return result


def plot(output, comparison):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = comparison['runs']
    labels = ['Order {}'.format(row['parameters']['order']) for row in rows]
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.2), constrained_layout=True)
    series = [([max(s['volume_vs_surface_relative_l2'] for s in row['checks']['surface_checks']
                    if s['order'] == max(t['order'] for t in row['checks']['surface_checks']))*100 for row in rows],
                'Volume vs closed surface'),
              ([row['checks']['two_surface_relative_l2']*100 for row in rows], 'Two enclosing surfaces'),
              ([row['checks']['volume_quadrature_relative_l2']*100 for row in rows], 'Volume quadrature change')]
    for ax, (values, title) in zip(axes, series):
        ax.bar(labels, values, color=['#547e9a', '#238c72'], width=.55)
        for i, value in enumerate(values):
            ax.text(i, value, '{:.3g}%'.format(value), ha='center', va='bottom')
        ax.set(title=title, ylabel='Relative complex-field L2 difference (%)', ylim=(0, max(values)*1.3))
        ax.grid(axis='y', alpha=.2)
        ax.set_axisbelow(True)
    fig.suptitle('Finite 3D grating: matched mesh controls | not yet converged', fontsize=12)
    fig.savefig(output/'convergence_diagnostics.png', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, required=True)
    parser.add_argument('--refined', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    result = compare(args.base, args.refined)
    (args.output/'comparison.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    plot(args.output, result)
    print(json.dumps({k: v for k, v in result.items() if k != 'runs'}, indent=2))


if __name__ == '__main__':
    main()
