"""Plot absolute collection and adjacent-step changes for an X-cell sequence."""
import argparse
from pathlib import Path

import numpy as np

from rcwa_prototype.current_series import save_report


def plot(result, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cells = result['cells']
    x = [c['lx_um'] for c in cells]
    values = np.array([c['collected_cross_section_um2'] for c in cells])
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), layout='constrained')
    for ch, (color, marker) in enumerate(zip(('#2878A5', '#BD5E34', '#35854D'), ('o', 's', '^'))):
        axes[0].plot(x, values[:, ch], marker=marker, color=color,
                     linestyle='--' if ch == 1 else '-', label='CH{}'.format(ch+1))
        axes[1].plot(x[1:], [100*step['channels'][ch]['relative_collected_power_change']
                            for step in result['adjacent_steps']], marker=marker,
                     color=color, linestyle='--' if ch == 1 else '-', label='CH{}'.format(ch+1))
    axes[0].set(ylabel='Continuous collected cross section (um^2)', title='Absolute power, no peak normalization')
    axes[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    axes[1].set(ylabel='Change from previous cell (%)', title='Adjacent steps, same continuous angles')
    for sign in (-1, 1):
        axes[1].axhline(sign*100*result['tolerances']['power_relative'], color='0.5', lw=.8, ls=':')
    axes[1].axhline(0., color='0.5', lw=.8)
    for ax in axes:
        ax.set(xlabel='X cell length (um)', xticks=x)
        ax.grid(alpha=.2)
        ax.legend(frameon=False)
    fig.suptitle('Periodic-current X-cell sensitivity, NOT an isolated-defect certification')
    fig.savefig(output/'x_cell_sequence.png', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radiation', type=Path, action='append', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = save_report(args.radiation, args.output)
    plot(result, args.output)
    lines = ['# X-cell sequence: periodic current only', '',
        '| Lx (um) | Mx / My | CH1 (um^2) | CH2 (um^2) | CH3 (um^2) |',
        '|---|---|---:|---:|---:|']
    for cell in result['cells']:
        lines.append('| {:.8g} | {} / {} | {:.9g} | {:.9g} | {:.9g} |'.format(
            cell['lx_um'], cell['harmonics_x'], cell['harmonics_y'], *cell['collected_cross_section_um2']))
    lines += ['', '| Step | Channel | Power change | Complex-field change |', '|---|---|---:|---:|']
    for step in result['adjacent_steps']:
        for ch, row in enumerate(step['channels'], 1):
            lines.append('| {} to {} | CH{} | {:+.6f}% | {:.6f}% |'.format(
                step['first_index']+1, step['second_index']+1, ch,
                100*row['relative_collected_power_change'], 100*row['relative_weighted_complex_field_change']))
    lines += ['', '| Step | Internal E / contrast-current change | Contrast moment change |',
              '|---|---:|---:|']
    for step in result['adjacent_steps']:
        current = step['internal_source']
        moment = current['contrast_moment_relative_change']
        lines.append('| {} to {} | {:.6f}% | {} |'.format(step['first_index']+1,
            step['second_index']+1, 100*current['internal_e_relative_l2'],
            '{:.6f}%'.format(100*moment) if moment is not None else 'N/A'))
    lines += ['', 'Last two local X steps passed: {}.'.format(result['last_two_x_steps_passed']),
              'Isolated defect certified: False.', '']+result['limitations']
    (args.output/'SERIES.md').write_text('\n'.join(lines), encoding='utf-8')
    print('\n'.join(lines), flush=True)


if __name__ == '__main__':
    main()
