"""Reviewable plots and CSV of numeric, statistical and footprint sensitivity."""
import csv


def plot_convergence(output, report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    axes_data = report['axes']
    fig, axes = plt.subplots(len(axes_data), 4, figsize=(16, 3*len(axes_data)),
                             layout='constrained', squeeze=False)
    rows = []
    for row, (name, scan) in enumerate(axes_data.items()):
        for col, observable in enumerate(('hemisphere', 'ch1', 'ch2', 'ch3')):
            means, errors = [], []
            for case in scan['cases']:
                stat = case['statistics']['observables'][observable]
                means.append(stat['mean_um2'])
                errors.append(stat['ci95_um2'][1]-stat['mean_um2'] if stat['ci95_um2'] else 0.)
                rows.append([name, case['level'], observable, stat['samples'], stat['mean_um2'],
                    stat['standard_error_um2'], *(stat['ci95_um2'] or [None, None]),
                    case['maximum_dofs'], case['peak_accounted_gib'], case['mean_sample_elapsed_seconds']])
            ax = axes[row, col]
            ax.errorbar([case['level'] for case in scan['cases']], means, yerr=errors,
                        marker='o', capsize=3, color=('tab:green', 'tab:red', 'tab:blue', 'tab:purple')[col])
            ax.set(xlabel=name, ylabel='DSC integral (um2)', title=observable)
            ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
            ax.grid(alpha=.25)
            if scan.get('limits'):
                ax.text(.02, .97, 'Incomplete levels: '+', '.join(
                    '{:g}'.format(item['level']) for item in scan['limits']),
                    transform=ax.transAxes, va='top', fontsize=8, color='firebrick')
    fig.suptitle('Haze convergence | 95% intervals | same-surface / shared-parent pairing')
    fig.savefig(output/'haze_convergence.png', dpi=140)
    plt.close(fig)
    plot_paired_errors(output, report)
    with (output/'haze_convergence.csv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['axis', 'level', 'observable', 'samples', 'mean_um2', 'standard_error_um2',
                         'ci95_low_um2', 'ci95_high_um2', 'maximum_dofs', 'peak_accounted_gib',
                         'mean_sample_seconds'])
        writer.writerows(rows)


def plot_paired_errors(output, report):
    import matplotlib.pyplot as plt
    import numpy as np

    scans = list(report['axes'].items())
    fig, axes = plt.subplots(1, 4, figsize=(16, max(4, .65*len(scans)+2)),
                             layout='constrained')
    target = report['inputs'].get('tolerance', .05)*100
    for ax, observable in zip(axes, ('hemisphere', 'ch1', 'ch2', 'ch3')):
        ax.axvline(target, color='black', linestyle='--', linewidth=1, label='Tolerance')
        for row, (_, scan) in enumerate(scans):
            comparison = scan.get('comparison') or {}
            item = comparison.get('observables', {}).get(observable, {})
            bound = item.get('relative_error_upper_bound')
            change = item.get('relative_change')
            if bound is not None:
                ax.barh(row, 100*bound, height=.5,
                        color='tab:green' if item['target_met'] else 'tab:red', alpha=.7)
            if change is not None:
                ax.plot(100*change, row, 'ko', markersize=4)
            if bound is None:
                ax.text(.02, row, comparison.get('status', 'unavailable'),
                        transform=ax.get_yaxis_transform(), va='center', fontsize=8)
        ax.set(yticks=np.arange(len(scans)), yticklabels=[name for name, _ in scans],
               xlabel='Change / paired upper bound\n(%, pointwise 95%)', title=observable)
        ax.set_xlim(left=0)
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=.25)
    fig.suptitle('Last two tested levels only | dots: mean change | bars: paired error bound')
    fig.savefig(output/'haze_paired_errors.png', dpi=140)
    plt.close(fig)
