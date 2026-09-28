"""Compare two fixed-geometry near audits with explicit interface clearances."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from .near_metrics import comparison
from .parameters import Parameters


def load(directory):
    directory = Path(directory)
    report = json.loads((directory/'report.json').read_text(encoding='utf-8'))
    report['parameters'] = asdict(Parameters(**report['parameters']))
    with np.load(directory/'fields.npz') as data:
        return report, dict(data)


def compare(base, refined, output):
    ra, a = load(base)
    rb, b = load(refined)
    for key, value in ra['parameters'].items():
        if key not in ('particle_mesh_um', 'pml_mesh_scale', 'order', 'pml_order', 'pml_inner_fraction') and value != rb['parameters'][key]:
            raise ValueError('Unexpected parameter change: '+key)
    for key, value in ra['roughness'].items():
        if key != 'mesh_um' and value != rb['roughness'][key]:
            raise ValueError('Roughness geometry differs: '+key)
    for key in ('points', 'zoom_points', 'triangles', 'surface_distance', 'sphere_distance'):
        np.testing.assert_allclose(a[key], b[key], atol=1e-14, rtol=0)
    whole = comparison(b['points'], a['electric'], b['electric'], np.ones(len(b['points']), bool))
    bulk = comparison(b['points'], a['electric'], b['electric'],
                      (b['surface_distance'] > .01) & (b['sphere_distance'] > .005))
    gap = []
    for clearance in (1, 2):
        mask = np.minimum(b['zoom_surface_distance'], b['zoom_sphere_distance']) >= clearance/1000
        gap.append(dict(clearance_nm=clearance, **comparison(b['zoom_points'], a['zoom_e'], b['zoom_e'], mask)))
    result = dict(base=ra, refined=rb, entire_yz_slice=whole, bulk=bulk, gap=gap,
        changed_parameters={key: [value, rb['parameters'][key]] for key, value in ra['parameters'].items()
                            if value != rb['parameters'][key]},
        changed_roughness={key: [value, rb['roughness'][key]] for key, value in ra['roughness'].items()
                           if value != rb['roughness'][key]},
        conclusion='Two discretizations, not an error bound or proof of asymptotic convergence; see changed_parameters. No smoothing. Whole-slice peaks include interfaces; gap metrics use fixed clearances.',
        gap_under_one_percent=all(c['sampled_peak_relative_change'] < .01 and c['complex_e_relative_l2'] < .01 for c in gap))
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output/'comparison.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    plot(result, a, b, output)
    print(json.dumps({k: v for k, v in result.items() if k not in ('base', 'refined')}, indent=2))


def plot(report, a, b, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    colors = ('#60757c', '#087f76')
    labels = []
    for key in ('base', 'refined'):
        p = report[key]['parameters']
        labels.append('p{}, PML{} (inner {:g}%), h x{:g}'.format(
            p['order'], p['pml_order'] or p['order'],
            100*p['pml_inner_fraction'] if p['pml_order'] else 100,
            report[key]['local_mesh_scale']))
    for data, label, color in zip((report['base'], report['refined']), labels, colors):
        rows = data['nested_slice_peaks']
        axes[0, 0].plot([r['samples'] for r in rows], [r['peak']['intensity'] for r in rows],
                        'o-', color=color, label=label)
        for clearance, style in ((1, 'o-'), (2, 's--')):
            rows = [r for r in data['gap_nested_peaks'] if r['clearance_nm'] == clearance]
            axes[0, 1].plot([r['samples'] for r in rows], [r['peak']['intensity'] for r in rows],
                            style, color=color, label=label+', clearance '+str(clearance)+' nm')
    for data, label, color in zip((a, b), labels, colors):
        intensity = np.sum(abs(data['zoom_e'])**2, axis=-1).reshape(len(data['zoom_z']), len(data['zoom_y']))
        axes[1, 0].plot(data['zoom_z']*1000, intensity[:, len(data['zoom_y'])//2], color=color, label=label)
    labels, rows = ['Whole slice', 'Bulk', 'Gap >=1 nm', 'Gap >=2 nm'], [report['entire_yz_slice'], report['bulk']]+report['gap']
    for offset, key, name, color in ((-.18, 'complex_e_relative_l2', 'Complex E L2', '#087f76'),
                                    (.18, 'sampled_peak_relative_change', 'Sampled peak', '#b35c37')):
        axes[1, 1].bar(np.arange(4)+offset, [100*r[key] for r in rows], width=.36, label=name, color=color)
    axes[0, 0].set(title='Whole YZ slice: sample-dependent peak', xlabel='Samples per axis', ylabel='Peak |E/E_inc|^2')
    axes[0, 1].set(title='Gap peak with explicit interface clearance', xlabel='Local samples per axis', ylabel='Peak |E/E_inc|^2')
    axes[1, 0].set(title='Fixed centerline x=y=0 inside gap', xlabel='Z (nm)', ylabel='|E/E_inc|^2')
    axes[1, 1].set(title='Discretization change: distinct error metrics', xticks=range(4), xticklabels=labels, ylabel='Relative change (%)')
    for ax in axes.flat:
        ax.legend(fontsize=8)
        ax.grid(alpha=.2)
        ax.set_axisbelow(True)
    fig.suptitle('Rough-particle near-field convergence: no smoothing or fitted normalization')
    fig.savefig(output/'near_convergence.png', dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', default='build/near_convergence/h70p2')
    parser.add_argument('--refined', default='build/near_convergence/base')
    parser.add_argument('--output', default='build/near_convergence/comparison')
    args = parser.parse_args()
    compare(args.base, args.refined, args.output)


if __name__ == '__main__':
    main()
