"""Report fresh Mie/FEM comparisons without treating literature-input runs as measurements."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import miepython
import numpy as np

from fullwave_prototype.parameters import Parameters
from fullwave_prototype.planar import background_vectors


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def mie_dsc(p, directions):
    down, _, _, ein, _, _ = background_vectors(p)
    frame = np.column_stack((ein, np.cross(down, ein), down))
    local = directions @ frame
    phi = np.arctan2(local[:, 1], local[:, 0])
    s1, s2 = miepython.S1_S2(p.particle_n, p.k0*p.diameter_nm/2000,
        np.clip(local[:, 2], -1, 1), norm='wiscombe')
    return (abs(s2*np.cos(phi))**2+abs(s1*np.sin(phi))**2)/p.k0**2


def angular_comparison(folder):
    report = read(folder/'report.json')
    p = Parameters(**report['parameters'])
    with np.load(folder/'fields.npz', allow_pickle=False) as raw:
        directions = raw['directions']
        actual = np.sum(abs(raw['amplitude'])**2, axis=-1)
    expected = mie_dsc(p, directions)
    n = int(np.sqrt(len(directions)/2))
    if 2*n*n != len(directions):
        raise ValueError('Expected the tensor sphere_rule quadrature')
    _, w = np.polynomial.legendre.leggauss(n)
    weights = np.repeat(w*np.pi/n, 2*n)
    reference_integral = float(weights@expected)
    if not np.isclose(reference_integral, report['exact']['csca_um2'], rtol=1e-10, atol=0):
        raise ValueError('Angular Mie normalization disagrees with its total cross section')
    np.savez_compressed(folder/'angular_comparison.npz', directions=directions,
        weights_sr=weights, fem_dsc_um2_sr=actual, mie_dsc_um2_sr=expected)
    return dict(relative_weighted_dsc_l2=float(np.sqrt(weights@(actual-expected)**2)/
        np.sqrt(weights@expected**2)), reference_integral_um2=reference_integral)


def figures(root, cases, literature):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), layout='constrained')
    for c in literature['cases']:
        with np.load(root/(c['name']+'.npz'), allow_pickle=False) as raw:
            t, y = raw['theta_signed_deg'], raw['dsc_l16_um2_sr']
        wavelength = c['simulation']['wavelength_nm']
        diameter = c['simulation']['particle_diameter_nm']
        if wavelength == 442:
            axes[0].semilogy(t, y, label='d = {:g} nm'.format(diameter))
        if diameter == 155:
            axes[1].semilogy(t, y, label='lambda = {:g} nm'.format(wavelength))
    for ax in axes:
        ax.set(xlim=(-90, 90), xlabel='Signed viewing angle (deg)', ylabel='DSC (um^2/sr)')
        ax.grid(alpha=.25)
        ax.legend()
    axes[0].set_title('442 nm | PSL on 1.5 nm SiO2 / Si')
    axes[1].set_title('155 nm PSL | coherent equal s/p')
    fig.suptitle('Recomputed paper inputs, NOT an overlay of experimental measurements')
    fig.savefig(root/'paper_inputs_recomputed.png', dpi=190)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), layout='constrained')
    for case in cases:
        report = read(root/case['folder']/'report.json')
        shell = report['near_shells']
        radius = [s['radius_nm'] for s in shell]
        label = '{} / h={:g} um'.format(case['folder'], report['parameters']['air_mesh_um'])
        axes[0].plot(radius, [100*s['electric_relative_l2'] for s in shell], 'o-', label=label)
        axes[1].plot(radius, [100*s['magnetic_relative_l2'] for s in shell], 'o-', label=label)
    for ax, name in zip(axes, ('Electric field E', 'Magnetic field Z0 H')):
        ax.set(xlabel='Probe-shell radius from sphere center (nm)', ylabel='Relative L2 error (%)', title=name)
        ax.grid(alpha=.25)
        ax.legend(fontsize=8)
    fig.suptitle('100 nm sphere in air, 266 nm | FEM versus independent Mie near fields')
    fig.savefig(root/'mie_fem_near_error.png', dpi=190)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    mie, nist = read(root/'mie_reference.json'), read(root/'nist_mie.json')
    literature = read(root/'literature_inputs_and_runs.json')
    exact = mie['direct_bessel']['csca_um2']
    with np.load(root/'nist_bv_air_raw.npz', allow_pickle=False) as raw:
        expected = mie_dsc(Parameters(**mie['parameters']), raw['directions'])
        actual, weights = raw['p_dsc_um2_sr'], raw['weights_sr']
        bv_air = dict(relative_weighted_dsc_l2=float(np.sqrt(weights@(actual-expected)**2)/
            np.sqrt(weights@expected**2)), nist_upper_dsc_um2=float(weights@actual),
            mie_upper_dsc_um2=float(weights@expected))
    cases, resources = [], []
    for path in sorted(root.glob('*.resources.json')):
        guard = read(path)
        resources.append(dict(folder=path.name.replace('.resources.json', ''),
            status=guard['status'], peak_accounted_gib=guard['peak_accounted_gib'],
            surviving_owned_pids=guard['surviving_owned_pids']))
        if guard['status'] != 'completed':
            continue
        folder = root/resources[-1]['folder']
        report = read(folder/'report.json')
        keys = ('wavelength_nm', 'diameter_nm', 'particle_n', 'particle_k',
            'substrate_n', 'substrate_k', 'gap_nm', 'incidence_deg', 'polarization',
            'beam_profile', 'beam_waist_um', 'defect_shape')
        if any(report['parameters'][key] != mie['parameters'][key] for key in keys):
            raise ValueError('FEM and Mie physical inputs differ: '+str(folder))
        angular = angular_comparison(folder)
        cases.append(dict(folder=folder.name, csca_um2=report['csca_um2'],
            signed_csca_error_percent=100*(report['csca_um2']/exact-1),
            extinction_error_percent=100*report['extinction_relative_error'],
            energy_balance_percent=100*report['energy_balance_relative'],
            angular_dsc_relative_l2_percent=100*angular['relative_weighted_dsc_l2'],
            max_near_e_l2_percent=100*max(s['electric_relative_l2'] for s in report['near_shells']),
            max_near_h_l2_percent=100*max(s['magnetic_relative_l2'] for s in report['near_shells']),
            dofs=report['solver']['dofs']))
    result = dict(mie_cross_section_um2=exact,
        bv_air_same_angle_polarization=bv_air,
        nist_mie_relative_error=abs(nist['csca_um2']/exact-1),
        independent_mie_relative_error=mie['scattering_relative_error'],
        fresh_fem_cases=cases, resources=resources,
        literature_measured_error=None, literature_status='matched_input_runs_only_not_digitized',
        caveats=['Homogeneous-air Mie validation is not a wafer-contact validation.',
            'No magnitude, position or global phase fitted in Mie/FEM comparison.',
            'No measured-data accuracy claim from multipole convergence.',
            'Legacy n=1.59,k=0 at 266 nm is an ideal control, not traced actual PSL dispersion.'])
    hashes = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob('*')) if p.is_file() and p.suffix in ('.npz', '.json')
        and p.name not in ('reference_results.json', 'data_manifest.json')}
    (root/'data_manifest.json').write_text(json.dumps(hashes, indent=2), encoding='utf-8')
    (root/'reference_results.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    figures(root, cases, literature)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
