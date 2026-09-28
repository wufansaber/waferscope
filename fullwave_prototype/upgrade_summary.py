"""Summarize completed physical experiments without certifying other cases."""
import argparse
import json
from pathlib import Path

import numpy as np


def read_report(root, case):
    return json.loads((root/case/'report.json').read_text(encoding='utf-8'))


def comparison(baseline, experiments, flat_results=None):
    improved = experiments/'order3_condensed'
    flat_results = flat_results or improved
    result = dict(scope='Specific 355 nm / 100 nm benchmark cases only; not a global accuracy certificate.',
                  mie=[], flat=[], pml=[], alternatives=[])
    for loss in ('0.0', '0.1'):
        before = read_report(baseline, 'mie_'+loss)
        after = read_report(improved, 'mie_'+loss)
        entry = dict(particle_k=float(loss))
        for name, report in (('before', before), ('after', after)):
            entry[name] = {key: report[key] for key in (
                'scattering_relative_error', 'extinction_relative_error', 'energy_balance_relative')}
            entry[name]['max_near_far_flux_relative_error'] = max(
                shell['scattered_vs_far_relative'] for shell in report['flux_shells'])
            entry[name]['all_listed_errors_below_1_percent'] = all(
                np.isfinite(value) and value < .01 for value in entry[name].values())
        result['mie'].append(entry)
    result['incomplete_checks'] = []
    for polarization in ('p', 's'):
        if (flat_results/('flat_'+polarization)/'report.json').exists():
            result['flat'].append(dict(polarization=polarization,
                source=str((flat_results/('flat_'+polarization)).resolve()),
                before=read_report(baseline, 'flat_'+polarization)['channels'],
                after=read_report(flat_results, 'flat_'+polarization)['channels']))
        else:
            result['incomplete_checks'].append('Order-3 flat-wafer '+polarization+' image validation')
    for name, root in (('before', baseline), ('after', improved)):
        cases = read_report(root, 'rough_pml')['cases']
        result['pml'].append(dict(configuration=name,
            max_power_change_vs_strength2=max(c['power_change_vs_strength2'] for c in cases),
            max_near_probe_change_vs_strength2=max(c['near_probe_relative_l2_vs_strength2'] for c in cases),
            cases=cases))
    for name in ('graded2', 'refined_pml', 'order3_condensed'):
        report = read_report(experiments/name, 'mie_0.0')
        result['alternatives'].append(dict(configuration=name,
            scattering_relative_error=report['scattering_relative_error'],
            energy_balance_relative=report['energy_balance_relative'], solver=report['solver']))
    return result


def plot_comparison(report, target):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    colors = ('#667781', '#087f77')
    labels = ('Previous order 2', 'Refined order 3')
    for index, key in enumerate(('before', 'after')):
        x = np.arange(2)+(index-.5)*.34
        axes[0, 0].bar(x, [100*c[key]['scattering_relative_error'] for c in report['mie']],
                       width=.34, color=colors[index], label=labels[index])
        axes[0, 1].bar(x, [100*c[key]['energy_balance_relative'] for c in report['mie']],
                       width=.34, color=colors[index])
        cases = report['pml'][index]
        axes[1, 0].bar(x, [100*cases['max_power_change_vs_strength2'],
                           100*cases['max_near_probe_change_vs_strength2']], width=.34, color=colors[index])
        complete_flat = report['flat'] if len(report['flat']) == 2 else []
        image_errors = [100*c['image_flux_relative_l2'] for pol in complete_flat for c in pol[key]]
        if image_errors:
            axes[1, 1].bar(np.arange(6)+(index-.5)*.34, image_errors, width=.34, color=colors[index])
    axes[0, 0].set(title='Sphere scattering vs Mie', xticks=[0, 1], xticklabels=['k=0', 'k=0.1'])
    axes[0, 0].legend()
    axes[0, 1].set(title='Sphere energy-balance discrepancy', xticks=[0, 1], xticklabels=['k=0', 'k=0.1'])
    axes[1, 0].set(title='Rough wafer: sensitivity to PML strength', xticks=[0, 1], xticklabels=['Hemisphere power', 'Near-field probes'])
    axes[1, 1].set(title='Smooth wafer images vs BV', xticks=np.arange(6),
                   xticklabels=['p +45', 'p -45', 'p 0', 's +45', 's -45', 's 0'])
    if len(report['flat']) != 2:
        axes[1, 1].text(.5, .5, 'Order-3 checks not completed\nLarge solves stopped after resource pressure\nNo new image-accuracy claim',
                       ha='center', va='center', transform=axes[1, 1].transAxes)
        axes[1, 1].set_xticks([])
    for ax in axes.flat:
        ax.set_ylabel('Relative difference (%)')
        ax.axhline(1, color='#b85445', linestyle='--', linewidth=1)
        ax.grid(axis='y', alpha=.2)
        ax.set_axisbelow(True)
    fig.suptitle('Measured accuracy checks, not a universal 1% certificate', fontsize=15)
    fig.savefig(target, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', default='build/physical_audit')
    parser.add_argument('--experiments', default='build/boundary_upgrade')
    parser.add_argument('--flat-results')
    parser.add_argument('--output')
    args = parser.parse_args()
    root = Path(args.experiments)
    report = comparison(Path(args.baseline), root, Path(args.flat_results) if args.flat_results else None)
    output = Path(args.output) if args.output else root
    output.mkdir(parents=True, exist_ok=True)
    (output/'comparison.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    plot_comparison(report, output/'accuracy_comparison.png')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
