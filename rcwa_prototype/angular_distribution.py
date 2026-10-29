"""Continuous upper-hemisphere radiation of a frozen cell contrast current."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from .current_radiation import load_source, matching_scenes
from .contrast_current import modal_reconstruction, project_ray


def hemisphere_grid(polar_samples, azimuthal_samples):
    for value, lower, upper in ((polar_samples, 4, 256), (azimuthal_samples, 8, 1024)):
        if isinstance(value, bool) or not isinstance(value, int) or not lower <= value <= upper:
            raise ValueError('Integer polar [4,256] and azimuthal [8,1024] samples required')
    nodes, weights = np.polynomial.legendre.leggauss(polar_samples)
    mu, wm = (nodes[::-1]+1)/2, weights[::-1]/2
    phi = (np.arange(azimuthal_samples)+.5)*2*np.pi/azimuthal_samples
    xx, pp = np.meshgrid(mu, phi, indexing='ij')
    rays = np.stack((np.sqrt(1-xx*xx)*np.cos(pp), np.sqrt(1-xx*xx)*np.sin(pp), xx), axis=-1)
    return dict(theta_deg=np.rad2deg(np.arccos(mu)), phi_deg=np.rad2deg(phi),
        directions=rays, solid_angle_weights=np.broadcast_to(wm[:, None]*2*np.pi/azimuthal_samples,
                                                             xx.shape).copy())


def evaluate(source, scene, polar_samples, azimuthal_samples, reference_halfwidth):
    arrays = hemisphere_grid(polar_samples, azimuthal_samples)
    if isinstance(reference_halfwidth, bool) or not isinstance(reference_halfwidth, int) or not 1 <= reference_halfwidth <= 128:
        raise ValueError('Reciprocal-background halfwidth must be an integer in [1,128]')
    fields = []
    for i, ring in enumerate(arrays['directions']):
        fields.append([project_ray(source, scene, ray, reference_halfwidth=reference_halfwidth) for ray in ring])
        print('CONTINUOUS_RING {}/{}'.format(i+1, polar_samples), flush=True)
    arrays['far_electric_um'] = np.asarray(fields)
    arrays['dsc_um2_per_sr'] = np.sum(abs(arrays['far_electric_um'])**2, axis=-1)
    if not all(np.isfinite(a).all() for a in arrays.values()):
        raise ValueError('Nonfinite continuous radiation; no completed result emitted')
    return arrays


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(source_folder, reference_folder, output, polar_samples=24, azimuthal_samples=72, reference_halfwidth=40):
    hemisphere_grid(polar_samples, azimuthal_samples)
    if isinstance(reference_halfwidth, bool) or not isinstance(reference_halfwidth, int) or not 1 <= reference_halfwidth <= 128:
        raise ValueError('Reciprocal-background halfwidth must be an integer in [1,128]')
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Use an empty output folder; existing results are never overwritten')
    source, scene, info, provenance = load_source(source_folder)
    reference, ref_scene, ref_info, ref_provenance = load_source(reference_folder)
    if not matching_scenes(scene, ref_scene) or not info.get('particle') or ref_info.get('particle'):
        raise ValueError('Matching particle/reference scenes required')
    if 'source_electric_fourier' not in source or not info.get('current_z_order'):
        raise ValueError('Select a contrast-current export, not an S4 modes-only scene')
    if scene.kind == 'grating' and scene.grating_method != 'analytic_li':
        raise ValueError('Analytic Li background required')
    for key in ('backend', 'backend_version', 'coordinates'):
        if key not in info or info[key] != ref_info.get(key):
            raise ValueError('Matching solver metadata required: '+key)
    if not np.array_equal(source['orders'], reference['orders']):
        raise ValueError('Matching source/reference modal orders required')
    rebuilt, gate = modal_reconstruction(source, scene, reference)
    errors = (gate.get('relative_power_weighted_complex_field_change'), gate.get('relative_total_power_difference'))
    if (not gate.get('comparable') or any(v is None or not np.isfinite(v) for v in errors)
            or errors[0] > 1e-7 or abs(errors[1]) > 2e-7):
        raise RuntimeError('Same-basis absolute complex-field gate failed')
    output.mkdir(parents=True, exist_ok=True)
    # Carry the source snapshots so exported bundles remain reproducible after moving them.
    import shutil
    for label, folder in (('source', source_folder), ('reference', reference_folder)):
        (output/label).mkdir()
        for filename in ('report.json', 'modes.npz'):
            shutil.copy2(Path(folder)/filename, output/label/filename)
    np.savez_compressed(output/'modal_gate.npz', reconstructed_electric=rebuilt,
        direct_electric=source['electric_reflected_z0']-reference['electric_reflected_z0'],
        directions=source['directions'], orders=source['orders'])
    arrays = evaluate(source, scene, polar_samples, azimuthal_samples, reference_halfwidth)
    np.savez_compressed(output/'angular_distribution.npz', **arrays)
    columns = np.column_stack((np.repeat(arrays['theta_deg'], azimuthal_samples),
        np.tile(arrays['phi_deg'], polar_samples), arrays['directions'].reshape(-1, 3),
        arrays['solid_angle_weights'].ravel(), arrays['far_electric_um'].real.reshape(-1, 3),
        arrays['far_electric_um'].imag.reshape(-1, 3), arrays['dsc_um2_per_sr'].ravel()))
    np.savetxt(output/'angular_distribution.csv', columns, delimiter=',',
        header='theta_deg,phi_deg,nx,ny,nz,solid_angle_sr,Fx_real_um,Fy_real_um,Fz_real_um,Fx_imag_um,Fy_imag_um,Fz_imag_um,DSC_um2_per_sr', comments='')
    report = dict(format='waferscope-continuous-current-radiation', version=1,
        status='FROZEN_PERIODIC_CURRENT_NOT_ISOLATED_DEFECT', scene=asdict(scene),
        source=provenance, reference=ref_provenance, same_basis_gate=gate,
        sampling=dict(polar=polar_samples, azimuthal=azimuthal_samples, reference_halfwidth=reference_halfwidth),
        hemisphere_cross_section_um2=float(np.sum(arrays['solid_angle_weights']*arrays['dsc_um2_per_sr'])),
        angular_convergence_checked=False, isolated_defect_certified=False,
        meaning='Continuous radiation from one cell particle contrast current in the bare wafer/grating background.',
        limitations=['Current includes periodic-copy interactions; not an isolated-defect solution.',
            'Not the total sample/reference angular intensity and not their intensity subtraction.',
            'Fourier cutoff, geometry slicing, reciprocal background and angular sampling need convergence checks.',
            'No finite-beam, finite-grating-edge or camera-DN prediction.'])
    from matplotlib.figure import Figure
    for three_d, filename in ((False, 'continuous_angles.png'), (True, 'continuous_3d.png')):
        fig = Figure(figsize=(11, 5))
        plot_distribution(fig, arrays, three_dimensional=three_d)
        fig.savefig(output/filename, dpi=160)
    report['files_sha256'] = {str(p.relative_to(output)).replace('\\', '/'): digest(p)
        for p in sorted(output.rglob('*')) if p.is_file()}
    # Completion marker is last; failed/cancelled jobs cannot be loaded as complete.
    (output/'continuous_report.json').write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    return report


def load_distribution(folder):
    folder = Path(folder).resolve()
    report = json.loads((folder/'continuous_report.json').read_text(encoding='utf-8'))
    if (report.get('format') != 'waferscope-continuous-current-radiation' or report.get('version') != 1
            or report.get('status') != 'FROZEN_PERIODIC_CURRENT_NOT_ISOLATED_DEFECT'):
        raise ValueError('Select a completed continuous radiation result')
    required = {'angular_distribution.npz', 'angular_distribution.csv', 'modal_gate.npz',
                'source/modes.npz', 'source/report.json', 'reference/modes.npz', 'reference/report.json'}
    if not required.issubset(report.get('files_sha256', {})):
        raise ValueError('Incomplete continuous radiation manifest')
    for relative, expected in report['files_sha256'].items():
        path = (folder/relative).resolve()
        try:
            path.relative_to(folder)
        except ValueError:
            raise ValueError('Invalid bundle path')
        if digest(path) != expected:
            raise ValueError('Continuous radiation hash check failed: '+relative)
    with np.load(folder/'angular_distribution.npz', allow_pickle=False) as raw:
        arrays = dict(raw)
    return report, arrays


def plot_distribution(fig, arrays, logarithmic=True, three_dimensional=False):
    from grating_scene.plotting import positive_norm
    fig.clear()
    dsc = arrays['dsc_um2_per_sr']
    norm = positive_norm([dsc], logarithmic)
    if three_dimensional:
        ax = fig.add_subplot(111, projection='3d')
        rays = arrays['directions']
        # Close the azimuth seam explicitly; radius encodes direction, color encodes DSC.
        rays = np.concatenate((rays, rays[:, :1]), axis=1)
        colors = np.concatenate((dsc, dsc[:, :1]), axis=1)
        from matplotlib import cm
        ax.plot_surface(*rays.transpose(2, 0, 1), facecolors=cm.viridis(norm(colors)),
                        rstride=1, cstride=1, shade=False)
        im = cm.ScalarMappable(norm=norm, cmap='viridis')
        ax.set(xlabel='kx/k', ylabel='ky/k', zlabel='kz/k', title='Continuous upper hemisphere')
        ax.set_box_aspect((2, 2, 1))
    else:
        ax = fig.add_subplot(121)
        im = ax.pcolormesh(arrays['phi_deg'], arrays['theta_deg'], dsc,
                           norm=norm, cmap='viridis', shading='nearest', rasterized=True)
        ax.set(xlabel='Azimuth phi (deg), +X toward +Y', ylabel='Polar theta from +Z (deg)',
               xlim=(0, 360), ylim=(0, 90), title='Continuous DSC (absolute scale)')
        cut = fig.add_subplot(122)
        for target in (0, 90, 180, 270):
            index = np.argmin(abs(arrays['phi_deg']-target))
            cut.plot(arrays['theta_deg'], dsc[:, index], label='phi={:.2f} deg'.format(arrays['phi_deg'][index]))
        if logarithmic and np.any(dsc > 0):
            cut.set_yscale('log')
        cut.set(xlabel='Polar theta (deg)', ylabel='DSC (um^2/sr)', title='Nearest sampled azimuth cuts')
        cut.legend(fontsize=8)
        cut.grid(alpha=.2)
    fig.colorbar(im, ax=ax, shrink=.75, label='DSC (um^2/sr)')
    fig.suptitle('One-cell contrast-current radiation | periodic interactions retained; isolation unchecked', fontsize=10)
    fig.tight_layout()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('source', 'reference', 'output'):
        parser.add_argument('--'+key, required=True, type=Path)
    parser.add_argument('--polar-samples', type=int, default=24)
    parser.add_argument('--azimuthal-samples', type=int, default=72)
    parser.add_argument('--reference-halfwidth', type=int, default=40)
    parser.add_argument('--max-memory-gib', type=float, default=8.)
    parser.add_argument('--timeout', type=float, default=3600.)
    parser.add_argument('--cancel-file', type=Path)
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    if args.worker:
        run(args.source, args.reference, args.output, args.polar_samples, args.azimuthal_samples, args.reference_halfwidth)
    else:
        from fullwave_prototype.guarded_run import run_guarded
        result = run_guarded([sys.executable, '-m', 'rcwa_prototype.angular_distribution', '--worker']+sys.argv[1:],
            args.output.with_suffix('.resources.json'), maximum=args.max_memory_gib, reserve=3.,
            timeout=args.timeout, cancel_file=args.cancel_file)
        if result['status'] != 'completed':
            raise SystemExit(1)


if __name__ == '__main__':
    main()
