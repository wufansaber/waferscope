"""Report fixed-cell single-axis cutoff sensitivity from original complex data."""
import argparse
import json
from pathlib import Path

import numpy as np

from rcwa_prototype.cutoff_check import compare_cutoff


def plot(result, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.3), layout='constrained')
    x = np.arange(len(result['channels']))
    colors = ['#2878A5', '#BD5E34', '#35854D'][:len(x)]
    for shift, key, filled in ((-.18, 'baseline', False), (.18, 'refined', True)):
        case = result[key]
        label = 'Mx={}, My={}'.format(*case['harmonics'])
        axes[0].bar(x+shift, case['cross_section_um2'], width=.34,
                    color=colors if filled else 'white', edgecolor=colors, label=label)
    axes[0].set(title='Unnormalized collection', ylabel='Collected cross section (um^2)')
    axes[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    axes[0].legend(frameon=False, fontsize=9)
    for ax, key, title, tolerance in (
            (axes[1], 'relative_collected_power_change', 'Power change (%)', 'power_relative'),
            (axes[2], 'relative_weighted_complex_field_change', 'Complex-field change (%)', 'complex_field_relative')):
        values = [100*row[key] if row[key] is not None else np.nan for row in result['channels']]
        bars = ax.bar(x, values, color=colors, width=.52)
        ax.bar_label(bars, fmt='%.3f', padding=3, fontsize=9)
        threshold = 100*result['tolerances'][tolerance]
        ax.axhline(threshold, color='0.5', ls=':', lw=1)
        if key == 'relative_collected_power_change':
            ax.axhline(-threshold, color='0.5', ls=':', lw=1)
        ax.axhline(0., color='0.5', lw=.6)
        ax.set(title=title)
        ax.margins(y=.25)
    for ax in axes:
        ax.set(xticks=x, xticklabels=['CH{}'.format(i+1) for i in x])
        ax.grid(axis='y', alpha=.2)
        ax.set_axisbelow(True)
    fig.suptitle('Fixed periodic cell: {} cutoff only, not isolated-defect certification'.format(
        result['refined_axis'].upper()))
    fig.savefig(output/'cutoff_sensitivity.png', dpi=180)
    plt.close(fig)


def save_report(baseline, refined, output):
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError('Use an empty output directory to preserve prior comparisons')
    result = compare_cutoff(baseline, refined)
    output.mkdir(parents=True, exist_ok=True)
    (output/'comparison.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    lines = ['# Fixed-cell Fourier-cutoff control', '',
        'Cell: {} um. Changed axis: {}.'.format(result['cell_um'], result['refined_axis']),
        'Baseline Mx/My: {}. Refined Mx/My: {}.'.format(
            result['baseline']['harmonics'], result['refined']['harmonics']), '',
        '| Channel | Baseline (um^2) | Refined (um^2) | Power change | Complex-field change |',
        '|---|---:|---:|---:|---:|']
    for i, row in enumerate(result['channels'], 1):
        lines.append('| CH{} | {:.9g} | {:.9g} | {} | {} |'.format(i,
            row['first_cross_section_um2'], row['second_cross_section_um2'],
            '{:+.6f}%'.format(100*row['relative_collected_power_change'])
            if row['relative_collected_power_change'] is not None else 'N/A',
            '{:.6f}%'.format(100*row['relative_weighted_complex_field_change'])
            if row['relative_weighted_complex_field_change'] is not None else 'N/A'))
    current = result['internal_source']
    background = result['background_control']
    lines += ['', 'Internal E relative L2: {:.6f}%.'.format(100*current['internal_e_relative_l2']),
        'Contrast moment relative change: {}.'.format('{:.6f}%'.format(100*current['contrast_moment_relative_change'])
            if current['contrast_moment_relative_change'] is not None else 'N/A'),
        'Local single step passed: {}.'.format(result['local_step_passed']),
        'Isolated defect certified: False.', '',
        'Bare-reference reflectance: {:.12g} to {:.12g}.'.format(
            background['baseline_reflectance'], background['refined_reflectance']),
        'Bare-reference complex-field control: {}.'.format(
            background.get('relative_power_weighted_complex_field_change', background.get('reason'))),
        background['norm_denominator'], '']+result['limitations']
    (output/'COMPARISON.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    plot(result, output)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--refined', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = save_report(args.baseline, args.refined, args.output)
    print(json.dumps(result, indent=2, allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
