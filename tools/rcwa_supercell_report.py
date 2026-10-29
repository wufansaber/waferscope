"""Audit periodic-cell size with fixed physical cutoffs and per-particle units."""
import argparse
import json
from pathlib import Path

from rcwa_prototype.comparison import compare_pair
from rcwa_prototype.geometry import Scene
from rcwa_prototype.supercell import compare_cells
from s4_prototype.solver import ROOT
from tools.rcwa_crosscheck_report import checked_load


def pair(paths):
    info, sample = checked_load(ROOT/paths['sample'])
    reference_info, reference = checked_load(ROOT/paths['reference'])
    if not info['particle'] or reference_info['particle'] or info['scene'] != reference_info['scene']:
        raise ValueError('Matched sample/reference geometry and particle flags required')
    if info['backend'] != reference_info['backend']:
        raise ValueError('Matched sample/reference backend required')
    for key in ('backend_version', 'polarization', 'factorization', 'executable_sha256',
                'epsilon_method', 'stripe_rule', 'particle_method'):
        if info.get(key) != reference_info.get(key):
            raise ValueError('Matched sample/reference numerical settings required: '+key)
    summary = compare_pair(reference, sample, Scene(**info['scene']))
    return dict(info=info, summary=summary)


def plot(output, result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    x = np.arange(len(result['channels']))
    fig, ax = plt.subplots(figsize=(8, 4.8), layout='constrained')
    for offset, prefix, color in ((-.17, 'coarse', '#007D8A'), (.17, 'fine', '#CE6939')):
        cell = result[prefix+'_cell_um']
        ax.bar(x+offset, [ch[prefix+'_periodic_cross_section_um2'] for ch in result['channels']],
               .34, label='Cell {:.4g} x {:.4g} um'.format(*cell), color=color)
    ax.set(xticks=x, xticklabels=['CH'+str(i+1) for i in x],
           ylabel='Periodic difference-field cross section (um^2)',
           title='Cell-size sensitivity, fixed physical Fourier cutoffs')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.legend(frameon=False)
    ax.grid(axis='y', alpha=.2)
    ax.set_axisbelow(True)
    fig.savefig(output/'supercell.png', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding='utf-8'))
    a, b = pair(config['coarse']), pair(config['fine'])
    result = dict(compare_cells(a, b), backend=a['info']['backend'], sources=config)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'supercell.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    plot(args.output, result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
