"""Separate numerical accuracy diagnostics from physical approximation errors."""
import argparse
import json
from pathlib import Path


def plot_local(output, report):
    from .imaging_plotting import _pyplot
    plt = _pyplot()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout='constrained')
    for order in sorted(set(c['fem_order'] for c in report['cases'])):
        cases = sorted((c for c in report['cases'] if c['fem_order'] == order), key=lambda c: c['height_scale'])
        axes[0].plot([c['input_stationary_rms_nm'] for c in cases],
                     [100*c['comparison']['complex_relative_l2'] for c in cases],
                     'o-', label='FEM order {}'.format(order))
    axes[0].axhline(1, color='gray', linestyle=':', label='1% reference')
    axes[0].set(xscale='log', yscale='log', xlabel='Stationary input RMS (nm)',
                ylabel='FEM vs SPM complex relative L2 (%)', title='Same surface shape, decreasing height')
    axes[0].legend()
    names = ['{:.3g}x / P{}'.format(c['height_scale'], c['fem_order']) for c in report['cases']]
    axes[1].plot(names, [100*c['source_quadrature_relative_l2_change'] for c in report['cases']],
                 'o-', color='seagreen', label='Source quadrature 5 vs 7')
    axes[1].plot(names, [max(100*c['batched_vs_native_background_relative_l2'], 1e-15)
                        for c in report['cases']], 's-', color='purple', label='Batch vs native background')
    axes[1].set(yscale='log', ylabel='Relative L2 (%)', title='Independent numerical checks')
    axes[1].legend()
    fig.suptitle('Haze only | fixed material inputs | 60 direction probes | NOT full-model certification')
    fig.savefig(Path(output)/'local_validation.png', dpi=150)
    plt.close(fig)


def plot_range(output, report):
    from .imaging_plotting import _pyplot
    plt = _pyplot()
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), layout='constrained')
    widths = [c['width_um'] for c in report['cases']]
    refined = report.get('independent_refinement')
    for i, color in enumerate(('teal', 'crimson', 'darkgoldenrod')):
        axes[0].plot(widths, [c['collected_spm_power_w'][i]*1000 for c in report['cases']],
                     'o-', color=color, label='Ch{}'.format(i+1))
        axes[1].plot(widths[1:], [100*c['from_previous_cap_complex_relative_l2'][i]
                                   for c in report['cases'][1:]], 'o-', color=color, label='Ch{}'.format(i+1))
        if refined:
            changes = refined['last_adjacent_complex_changes']
            checked = changes.get('192', changes.get(192))
            label = '192-order check' if i == 0 else None
            axes[0].scatter(widths[-1], refined['collected_spm_power_w_192'][i]*1000,
                            color=color, marker='*', s=90, label=label)
            axes[1].scatter(widths[-1], 100*checked[i], color=color, marker='*', s=90, label=label)
    axes[0].set(xlabel='Rough-surface width (um)', ylabel='Collected SPM power (mW)', title='Absolute power')
    axes[0].legend()
    axes[1].axhline(1, color='gray', linestyle=':', label='1% target')
    axes[1].set(xlabel='Rough-surface width (um)', ylabel='Adjacent weighted complex L2 (%)',
                title='Whole cap: phase-sensitive field change')
    axes[1].legend()
    axes[2].plot(widths, [100*c['input_power_coverage'] for c in report['cases']], 'o-', color='seagreen')
    axes[2].set(xlabel='Rough-surface width (um)', ylabel='Weighted incident coverage (%)',
                title='Coverage is NOT scattering accuracy')
    fig.suptitle('Same original synthesis parent | no new draw or periodic repetition | outer first-order SPM only')
    if refined:
        axes[0].set_title('Absolute power | original overview + refinement')
    fig.savefig(Path(output)/'range_validation.png', dpi=150)
    plt.close(fig)


def plot_mesh(output, report):
    from .imaging_plotting import _pyplot
    plt = _pyplot()
    cases = report['cases'][1:]
    names = [c['name'].replace('_', '\n') for c in cases]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), layout='constrained')
    for channel in range(3):
        for ax, key in zip(axes, ('complex_relative_l2', 'power_relative_changes')):
            ax.plot(names, [max(100*c[key][channel], 1e-8) for c in cases], 'o-', label='Ch{}'.format(channel+1))
    for ax, title in zip(axes, ('Complex field change', 'Collected power change')):
        ax.axhline(1, color='gray', linestyle=':', label='1% target')
        ax.set(yscale='log', ylabel='Relative change (%)', title=title)
        ax.legend()
    fig.suptitle('Same finite band-limited parent | {:.4g} um range | PSD tail not refined'.format(report['width_um']))
    fig.savefig(Path(output)/'surface_mesh.png', dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--local')
    parser.add_argument('--range')
    parser.add_argument('--mesh')
    args = parser.parse_args()
    for directory, name, plot in ((args.local, 'validation_report.json', plot_local),
                                   (args.range, 'range_report.json', plot_range),
                                   (args.mesh, 'surface_mesh_report.json', plot_mesh)):
        if directory:
            path = Path(directory)
            plot(path, json.loads((path/name).read_text(encoding='utf-8')))


if __name__ == '__main__':
    main()
