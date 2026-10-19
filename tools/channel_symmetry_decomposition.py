"""Expose particle-induced field, haze and signed interference on fixed scales."""
import argparse
import json
from pathlib import Path

import numpy as np

from .channel_symmetry_local import comparison


def run(source, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    names = ('change', 'without', 'interference', 'with')
    labels = ('Particle-induced intensity', 'Haze only', 'Signed interference', 'Total with particle')
    arrays, channels = {}, []
    with np.load(source/'image_fields.npz') as saved:
        for i in (1, 2):
            prefix = 'ch{}_'.format(i)
            axis = saved[prefix+'image_axis_um']
            center = np.argmin(abs(axis))
            selected = abs(axis) <= 12
            row = {}
            for name in names:
                field = saved[prefix+'image_'+name+'_w_um2']
                arrays[i, name] = field[np.ix_(selected, selected)]
                row[name] = dict(center_w_um2=float(field[center, center]),
                                 peak_in_24um_box_w_um2=float(arrays[i, name].max()))
            channels.append(row)
    report = dict(channels=channels,
        particle_induced_peak_comparison=comparison(
            channels[0]['change']['peak_in_24um_box_w_um2'],
            channels[1]['change']['peak_in_24um_box_w_um2']),
        scope='Particle-induced field is with-minus-without complex field, including coupling; '
              'it is not the total camera signal and not an isolated-particle Maxwell solution. '
              'The signed interference cannot be assigned a positive-intensity normalization.')
    (output/'peak_decomposition.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    fig, axes = plt.subplots(4, 2, figsize=(10, 15), constrained_layout=True)
    maximum = max(float(arrays[i, name].max()) for i in (1, 2) for name in names if name != 'interference')
    signed_max = max(float(abs(arrays[i, 'interference']).max()) for i in (1, 2))
    for row, (name, label) in enumerate(zip(names, labels)):
        for col in (0, 1):
            signed = name == 'interference'
            im = axes[row, col].pcolormesh(axis[selected], axis[selected], arrays[col+1, name],
                shading='auto', cmap='RdBu_r' if signed else 'inferno',
                vmin=-signed_max if signed else 0, vmax=signed_max if signed else maximum)
            axes[row, col].set(title='{} | CH{}'.format(label, col+1), aspect='equal',
                               xlabel='Image x (um)', ylabel='Image y (um)')
        fig.colorbar(im, ax=list(axes[row]), label='Signed W/um^2' if signed else 'W/um^2 | fixed absolute scale')
    fig.suptitle('No per-channel normalization | same coordinates | no detector noise')
    fig.savefig(output/'peak_decomposition.png', dpi=150)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(9, 8), constrained_layout=True)
    for row, (name, label) in enumerate((('change', 'Particle-induced field'), ('with', 'With haze + interference'))):
        for col in (0, 1):
            im = axes[row, col].pcolormesh(axis[selected], axis[selected], arrays[col+1, name],
                shading='auto', cmap='inferno', vmin=0, vmax=maximum)
            axes[row, col].set(title='{} | CH{}'.format(label, col+1), aspect='equal',
                               xlabel='Image x (um)', ylabel='Image y (um)')
    fig.colorbar(im, ax=list(axes.flat), label='W/um^2 | same absolute scale for all panels')
    fig.savefig(output/'peak_comparison.png', dpi=160)
    plt.close(fig)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('build/field_coverage_20261016/result'))
    parser.add_argument('--output', type=Path, default=Path('build/channel_symmetry_20261016'))
    args = parser.parse_args()
    run(args.source, args.output)
