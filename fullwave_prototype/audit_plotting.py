"""Render physical-validation evidence without presenting it as certification."""
import json
import os
from pathlib import Path

import numpy as np


def main():
    os.environ.setdefault('MPLCONFIGDIR', str(Path('build/fem_validation/.matplotlib').resolve()))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root = Path('build/physical_audit')

    def read(name):
        return json.loads((root/name/'report.json').read_text())

    cases = [read(name) for name in ('mie_0.0', 'mie_0.1', 'mie_order3')]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), layout='constrained')
    ax = axes[0, 0]
    x = np.arange(3)
    for offset, key, label in ((-.25, 'scattering_relative_error', 'Scattering vs Mie'),
                               (0, 'extinction_relative_error', 'Extinction vs Mie'),
                               (.25, 'energy_balance_relative', 'Extinction - scatter - absorption')):
        ax.bar(x+offset, [100*c[key] for c in cases], width=.25, label=label)
    ax.set_xticks(x, ['Order 2\nlossless', 'Order 2\nk=0.1', 'Order 3\nlossless'])
    ax.set(ylabel='Relative error (%)', title='Free sphere: independent Mie and energy balance')
    ax.legend(fontsize=8)
    ax.grid(axis='y', alpha=.2)

    with np.load(root/'airy.npz') as data:
        ax = axes[0, 1]
        ax.plot(data['axis_um'], data['airy'], label='Analytic Airy', color='black')
        ax.plot(data['axis_um'], data['calculated'], '--', label='Vector propagation, paraxial limit', color='#c45a25')
        ax.set(title='Independent ideal-optics check (no fitted amplitude)', xlabel='Image x (um)', ylabel='Unit-pupil axial flux')
        ax.legend(fontsize=8)

    pml = read('rough_pml')['cases']
    ax = axes[1, 0]
    ax.plot([c['strength'] for c in pml], [c['upper_dsc_um2'] for c in pml], 'o-', color='#237c70')
    ax.set(title='Rough particle: same mesh, varied PML strength', xlabel='PML strength', ylabel='Upper DSC (um^2)')
    ax.grid(alpha=.2)

    ax = axes[1, 1]
    colors = ('#24568c', '#c45a25', '#237c70')
    with np.load(root/'flat_p/images.npz') as data:
        xx = data['axis_um']
        roi = abs(xx) <= 20
        for i, color in enumerate(colors, 1):
            fem, bv = data['ch{}_fem'.format(i)], data['ch{}_bv'.format(i)]
            ax.plot(xx[roi], bv[len(xx)//2, roi], color=color, label='Ch{} BV'.format(i))
            ax.plot(xx[roi], fem[len(xx)//2, roi], '--', color=color, label='Ch{} FEM'.format(i))
    ax.set(title='Flat substrate: origin-aligned coherent image cuts', xlabel='Image x (um), y=0', ylabel='Image axial flux / incident irradiance')
    ax.legend(fontsize=8, ncol=2)
    fig.suptitle('PHYSICAL VALIDATION | Default order-2 FEM is not certified to 1% accuracy\n'
                 'Benchmarks constrain particular cases; they do not prove arbitrary rough-surface near-field accuracy', fontsize=12)
    fig.savefig(root/'validation_evidence.png', dpi=160)
    plt.close(fig)


if __name__ == '__main__':
    main()
