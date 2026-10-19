"""Compare grid masks and cone quadrature against saved FEM volume sources."""
import argparse
import json
from pathlib import Path
import time

import numpy as np

from fullwave_prototype.aperture_integrals import integrate_maps, integrate_sources
from fullwave_prototype.far_field import angular_grid, radiation_amplitude
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.run import hemisphere_integral


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    folder, output = Path(args.sample), Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    report = json.loads((folder/'sample_report.json').read_text(encoding='utf-8'))
    p = Parameters(**report['parameters'])
    apertures, tilts = (.5, .5, .5), (45, -45, 0)
    with np.load(folder/'haze_sources.npz', allow_pickle=False) as data:
        points, sources = data['points'], data['weights']
    start = time.perf_counter()
    direct = {}
    for order in (8, 16, 32):
        print('Direct FEM-source cap quadrature order '+str(order), flush=True)
        direct[str(order)] = integrate_sources(points, sources, p, apertures, tilts, order).tolist()
    reference = np.array(direct['32'])
    grids = {}
    for step in (10, 5, 2.5):
        print('Angular map step '+str(step), flush=True)
        theta, phi, rays = angular_grid(step)
        amplitude = radiation_amplitude(points, sources, rays.reshape(-1, 3), p)
        intensity = np.sum(abs(amplitude.reshape(rays.shape[:2]+(2,)))**2, axis=-1)
        modern = integrate_maps(intensity[None], theta, phi, apertures, tilts)[0]
        legacy = []
        for na, tilt in zip(apertures, tilts):
            axis = np.array((np.sin(np.deg2rad(tilt)), 0, np.cos(np.deg2rad(tilt))))
            mask = rays @ axis >= np.sqrt(1-na*na)
            legacy.append(hemisphere_integral(intensity*mask, theta))
        grids[str(step)] = dict(mask_um2=legacy, cone_um2=modern.tolist(),
            mask_relative_error=(np.abs(np.array(legacy)/reference-1)).tolist(),
            cone_relative_error=(np.abs(modern/reference-1)).tolist())
    result = dict(sample=str(folder.resolve()), parameters=report['parameters'],
        surface_sha256=report['surface_sha256'], numerical_apertures=list(apertures),
        tilt_degrees=list(tilts), direct_source_orders_um2=direct, angular_grids=grids,
        direct_16_to_32_relative_change=(np.abs(np.array(direct['16'])/reference-1)).tolist(),
        elapsed_seconds=time.perf_counter()-start,
        scope='One saved FEM field; independent angular evaluation, not independent Maxwell accuracy')
    (output/'aperture_audit.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), layout='constrained')
    for index, ax in enumerate(axes):
        for method, color in (('mask', 'tab:red'), ('cone', 'tab:green')):
            ax.plot([10, 5, 2.5], [grids[str(step)][method+'_um2'][index] for step in (10, 5, 2.5)],
                    'o-', color=color, label=method)
        ax.axhline(reference[index], color='black', linestyle='--', label='FEM source, order 32')
        ax.set(title='Ch'+str(index+1), xlabel='Angular step (deg)', ylabel='DSC integral (um2)')
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
        ax.grid(alpha=.25)
    axes[0].legend(fontsize=8)
    fig.suptitle('Objective solid-angle integration | same FEM volume sources')
    fig.savefig(output/'aperture_audit.png', dpi=150)
    plt.close(fig)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
