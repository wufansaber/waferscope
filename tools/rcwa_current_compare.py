"""Report periodic-current cell sensitivity after separating angular sampling."""
import argparse
import json
from pathlib import Path

import numpy as np

from rcwa_prototype.current_compare import compare_folders


def plot(output, result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout='constrained')
    x = np.arange(len(result['channels']))
    for offset, key, cell, color in ((-.17, 'first_cross_section_um2', result['first_lx_um'], '#4878A8'),
                                    (.17, 'second_cross_section_um2', result['second_lx_um'], '#CF8050')):
        axes[0].bar(x+offset, [r[key] for r in result['channels']], .34,
                    label='Lx={:.4g} um'.format(cell), color=color)
    axes[0].set(ylabel='Continuous collected cross section (um^2)',
                title='Separately solved periodic currents, NA=0.6')
    axes[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    for offset, key, label, color in ((-.17, 'relative_collected_power_change', 'Power change (signed)', '#CF8050'),
                                    (.17, 'relative_weighted_complex_field_change', 'Complex-field change', '#4878A8')):
        values = [100*r[key] if r[key] is not None else np.nan for r in result['channels']]
        axes[1].bar(x+offset, values, .34, label=label, color=color)
    axes[1].set(ylabel='Relative change (%)', title='Identical angles and absolute phase reference')
    axes[1].axhline(0., color='0.4', lw=.8)
    for ax in axes:
        ax.set(xticks=x, xticklabels=['CH{}'.format(i+1) for i in x])
        ax.grid(axis='y', alpha=.2)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=9)
    fig.suptitle('Cell sensitivity remains: not an isolated-defect certification')
    fig.savefig(output/'cell_sensitivity.png', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first', required=True, type=Path)
    parser.add_argument('--second', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        parser.error('Use an empty output directory to preserve previous comparisons')
    result = compare_folders(args.first, args.second)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'comparison.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    lines = ['# Continuous-angle cell sensitivity', '',
             '| Channel | First cross section (um^2) | Second (um^2) | Power change | Complex-field change |',
             '|---|---:|---:|---:|---:|']
    for i, row in enumerate(result['channels'], 1):
        power = row['relative_collected_power_change']
        field = row['relative_weighted_complex_field_change']
        lines.append('| CH{} | {:.9g} | {:.9g} | {} | {} |'.format(i,
            row['first_cross_section_um2'], row['second_cross_section_um2'],
            '{:+.6f}%'.format(100*power) if power is not None else 'N/A',
            '{:.6f}%'.format(100*field) if field is not None else 'N/A'))
    lines += ['', result['caveat'], '']
    (args.output/'COMPARISON.md').write_text('\n'.join(lines), encoding='utf-8')
    plot(args.output, result)
    print(json.dumps(result['channels'], indent=2))


if __name__ == '__main__':
    main()
