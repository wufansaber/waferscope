"""Compare two actual large-cross FEM previews without amplitude/phase fitting."""
import argparse
import json
from pathlib import Path

import numpy as np


def load(folder):
    report = json.loads((folder/'imaging_report.json').read_text(encoding='utf-8'))
    source = Path(report['source_directory'])
    scattering = json.loads((source/'report.json').read_text(encoding='utf-8'))
    with np.load(folder/'image_fields.npz') as saved:
        images = [saved['ch{}_image_with_w_um2'.format(i)] for i in (1, 2, 3)]
        axes = [saved['ch{}_image_axis_um'.format(i)] for i in (1, 2, 3)]
    return report, scattering, images, axes


def check(a, b, output):
    old, new = load(a), load(b)
    for key in ('wavelength_nm', 'defect_shape', 'defect_length_nm', 'defect_width_nm',
                'defect_height_nm', 'defect_rotation_deg', 'particle_n', 'particle_k',
                'substrate_n', 'substrate_k', 'incidence_deg', 'polarization', 'gap_nm', 'domain_radius_um'):
        if old[0]['scattering_parameters'][key] != new[0]['scattering_parameters'][key]:
            raise ValueError('Different physics: '+key)
    rows = []
    for i in range(3):
        np.testing.assert_allclose(old[3][i], new[3][i], atol=1e-12)
        selected = abs(new[3][i]) <= 32
        ia, ib = (case[2][i][np.ix_(selected, selected)] for case in (old, new))
        rows.append(dict(channel=i+1, image_relative_l2=float(np.linalg.norm(ia-ib)/np.linalg.norm(ib)),
            power_relative_change=float(new[0]['channels'][i]['with']['image_power_w']/
                                        old[0]['channels'][i]['with']['image_power_w']-1)))
    def symmetry(case):
        image, axis = case[2][2], case[3][2]
        selected = abs(axis) <= 32
        mirror = (-np.arange(len(axis))) % len(axis)
        x = image[np.ix_(selected, selected)]
        y = image[:, mirror][np.ix_(selected, selected)]
        power = [ch['with']['image_power_w'] for ch in case[0]['channels']]
        return dict(ch3_x_mirror_relative_l2=float(np.linalg.norm(x-y)/np.linalg.norm(x)),
                    side_power_relative_difference=float(2*(power[0]-power[1])/(power[0]+power[1])))
    result = dict(coarse=str(a.resolve()), refined=str(b.resolve()), channels=rows,
        coarse_symmetry=symmetry(old), refined_symmetry=symmetry(new),
        dofs=[case[1]['metrics']['with']['dofs'] for case in (old, new)],
        scope='Background mesh sensitivity only; local defect mesh fixed at 30 nm. '
              'No normalization, alignment, symmetrization or fitted scale. '
              'Not a full FEM/PML/local-mesh convergence proof.')
    output.mkdir(parents=True, exist_ok=True)
    (output/'mesh_sensitivity.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import PowerNorm
    maximum = max(image.max() for case in (old, new) for image in case[2])
    norm = PowerNorm(.5, vmin=0, vmax=maximum)
    fig, axes = plt.subplots(2, 3, figsize=(14, 9), constrained_layout=True)
    for row, case in enumerate((old, new)):
        for column, (image, axis) in enumerate(zip(case[2], case[3])):
            mask = abs(axis) <= 34
            im = axes[row, column].pcolormesh(axis[mask], axis[mask], image[np.ix_(mask, mask)],
                                             shading='auto', cmap='magma', norm=norm)
            axes[row, column].set(title='{} | CH{} | {} DOFs'.format(
                'Coarse background' if row == 0 else 'Finer background', column+1, result['dofs'][row]),
                xlabel='Image x (um)', ylabel='Image y (um)', xlim=(-32, 32), ylim=(-32, 32), aspect='equal')
    fig.colorbar(im, ax=list(axes.flat), label='W / um^2 | one absolute sqrt scale', shrink=.8)
    fig.suptitle('2 um CROSS | Actual FEM background-mesh sensitivity | no image symmetrization\n'
                 'Qualitative shape only: these preview results are not fully converged')
    fig.savefig(output/'mesh_sensitivity.png', dpi=140)
    plt.close(fig)
    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coarse', type=Path, default=Path('build/cross_2um_20261017/preview_scale3'))
    parser.add_argument('--refined', type=Path, default=Path('build/cross_2um_20261017/background_refined'))
    parser.add_argument('--output', type=Path, default=Path('build/cross_2um_20261017'))
    args = parser.parse_args()
    check(args.coarse, args.refined, args.output)


if __name__ == '__main__':
    main()
