"""Source-backed plots and a Chinese verification note for the sphere-size audit."""
import argparse
import json
from pathlib import Path

import numpy as np

from tools.particle_brightness_metrics import local_exponents


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def require_completed(folder):
    resource = folder.parent/(folder.name+'.resources.json')
    if not resource.is_file() or read(resource).get('status') != 'completed':
        raise ValueError('FEM comparison requires a completed guarded run: '+str(folder))


def series(rows, incidence, polarization, channel, metric):
    selected = sorted((r for r in rows if r['incidence_deg'] == incidence and
                      r['polarization'] == polarization and r['channel'] == channel),
                      key=lambda r: r['diameter_nm'])
    return np.array([r['diameter_nm'] for r in selected]), np.array([r[metric] for r in selected])


def side_symmetry(rows):
    result = []
    for incidence in (0., 75.):
        for polarization in ('p', 's'):
            record = dict(incidence_deg=incidence, polarization=polarization,
                          definition='2*abs(CH1-CH2)/(CH1+CH2); no per-channel normalization')
            for key in ('cap_dsc_um2', 'pupil_dsc_um2', 'expected_peak_signal_dn'):
                d1, a = series(rows, incidence, polarization, 1, key)
                d2, b = series(rows, incidence, polarization, 2, key)
                if not np.array_equal(d1, d2) or not len(a) or np.any(a+b <= 0):
                    raise ValueError('Matching positive side-channel series required')
                record['max_relative_'+key+'_asymmetry'] = float(np.max(2*abs(a-b)/(a+b)))
            result.append(record)
    return result


def mesh_convergence(base, names=None):
    from tools.particle_brightness_phase import constant_phase_relative_l2
    if names is None:
        names = ('fem_order2_r022_pml2', 'fem_order2_r022_pml2_m0875', 'fem_order2_r022_pml2_m075')
    result = []
    for old_name, new_name in zip(names[:-1], names[1:]):
        old_dir, new_dir = (base/name/'d100' for name in (old_name, new_name))
        if not all((path/'imaging_report.json').is_file() for path in (old_dir, new_dir)):
            continue
        old, new = (read(path/'imaging_report.json') for path in (old_dir, new_dir))
        stable = lambda report: {k: v for k, v in report['scattering_parameters'].items()
                                 if k not in ('air_mesh_um', 'substrate_mesh_um')}
        if stable(old) != stable(new) or old['settings'] != new['settings']:
            raise ValueError('Mesh convergence requires matching physical and imaging inputs')
        for folder, report in ((old_dir, old), (new_dir, new)):
            require_completed(folder)
            if len(report['channels']) != 3:
                raise ValueError('Mesh convergence requires all three channels')
        with np.load(old_dir/'image_fields.npz', allow_pickle=False) as a, \
                np.load(new_dir/'image_fields.npz', allow_pickle=False) as b:
            for index, (c, d) in enumerate(zip(old['channels'], new['channels']), 1):
                power = d['with']['image_power_w']/c['with']['image_power_w']-1
                field = constant_phase_relative_l2(b['ch{}_pupil_with'.format(index)],
                                                  a['ch{}_pupil_with'.format(index)])
                result.append(dict(previous=old_name, current=new_name, channel=index,
                    substrate_mesh_um=new['scattering_parameters']['substrate_mesh_um'],
                    pml_strength=new['scattering_parameters']['pml_strength'],
                    relative_power_change=power, pupil_relative_l2=field,
                    within_step_targets=bool(abs(power) < .01 and field < .02)))
    return result


def comparisons(base, nist, nist_name='nist'):
    convergence = []
    for path in sorted((base/nist_name/'convergence').glob('*/report.json')):
        test = read(path)
        d = test['simulation']['particle_diameter_nm']
        for row in test['rows']:
            original = next(r for r in nist['rows'] if r['diameter_nm'] == d and r['incidence_deg'] == row['incidence_deg'] and
                            r['polarization'] == row['polarization'] and r['channel'] == row['channel'])
            convergence.append(dict(case=path.parent.name, channel=row['channel'], polarization=row['polarization'],
                incidence_deg=row['incidence_deg'], diameter_nm=d, pupil_samples=row['pupil_samples'],
                lmax=row['lmax'], cap_order=row.get('cap_order', 24),
                cap_change=row['cap_dsc_um2']/original['cap_dsc_um2']-1,
                pupil_change=row['pupil_dsc_um2']/original['pupil_dsc_um2']-1,
                camera_peak_change=row['expected_peak_signal_dn']/original['expected_peak_signal_dn']-1))
    fem = []
    for path in sorted(base.glob('fem*/d*/imaging_report.json')):
        report = read(path)
        p = report['scattering_parameters']
        if (p['defect_shape'] != 'sphere' or p['incidence_deg'] not in (0, 75)
                or p['wavelength_nm'] != 266 or p['gap_nm'] != 0 or p['polarization'] != 'p'
                or p['beam_profile'] != 'plane' or report.get('roughness', {}).get('rms_nm', 0) != 0):
            raise ValueError('FEM comparison requires actual 266-nm contact spheres at a tested incidence')
        for key in ('particle_n', 'particle_k', 'substrate_n', 'substrate_k'):
            if p[key] != nist['controls']['parameters'][key]:
                raise ValueError('FEM and NIST material inputs differ: '+key)
        for key in ('numerical_apertures', 'tilt_degrees', 'magnification', 'transmission', 'irradiance_w_um2'):
            if not np.allclose(report['settings'][key], nist['controls']['settings'][key], rtol=1e-12, atol=0):
                raise ValueError('FEM and NIST imaging inputs differ: '+key)
        require_completed(path.parent)
        with np.load(path.parent/'image_fields.npz', allow_pickle=False) as raw:
            for index, ch in enumerate(report['channels'], 1):
                original = next(r for r in nist['rows'] if r['diameter_nm'] == p['diameter_nm'] and
                                r['incidence_deg'] == p['incidence_deg'] and r['polarization'] == 'p' and r['channel'] == index)
                peak = float(raw['ch{}_detector_mean_signal_e'.format(index)].max()/report['settings']['detector']['gain_e_dn'])
                fem.append(dict(diameter_nm=p['diameter_nm'], incidence_deg=p['incidence_deg'], channel=index, fem_peak_dn=peak,
                    nist_peak_dn=original['expected_peak_signal_dn'],
                    relative_peak_change=peak/original['expected_peak_signal_dn']-1,
                    relative_power_change=ch['with']['image_power_w']/original['image_power_w']-1,
                    order=p['order'], source=str(path),
                    numerics={key: p[key] for key in ('domain_radius_um', 'pml_order', 'pml_power',
                              'pml_mesh_scale', 'pml_strength', 'air_mesh_um', 'substrate_mesh_um',
                              'particle_mesh_um', 'integration_bonus')}))
    return convergence, fem


def figures(base, nist, nist_name='nist'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = nist['rows']
    colors = ('#2C718E', '#CF8050', '#3A9275')
    fig, all_axes = plt.subplots(2, 3, figsize=(13, 7.8), layout='constrained')
    for incidence, axes in zip((0., 75.), all_axes):
        for channel, color in enumerate(colors, 1):
            d, sigma = series(rows, incidence, 'p', channel, 'cap_dsc_um2')
            _, peak = series(rows, incidence, 'p', channel, 'expected_peak_signal_dn')
            axes[0].loglog(d, sigma, 'o-', label='CH{}'.format(channel), color=color)
            axes[1].loglog(d, peak, 'o-', label='CH{}'.format(channel), color=color)
            axes[2].plot(np.sqrt(d[:-1]*d[1:]), local_exponents(d, sigma), 'o-', label='CH{}'.format(channel), color=color)
        axes[2].axhline(6., color='0.4', ls=':', label='Rayleigh exponent 6')
        axes[0].set(ylabel='Collected cross section (um^2)', title='{:g} deg P: absolute collection'.format(incidence))
        axes[1].set(ylabel='Expected peak net signal (DN)', title='Pixel-integrated, black offset excluded')
        axes[2].set(ylabel='Local slope d log(P) / d log(d)', title='Interface-modified size law')
        for ax in axes:
            ax.set_xlabel('Sphere diameter (nm)')
            ax.grid(alpha=.2)
            ax.legend(frameon=False, fontsize=8)
    fig.suptitle('Contact spheres, gap=0 | 266 nm | NA=0.5 | 20x | 10 us | T=1e-6 | fixed incident irradiance')
    fig.savefig(base/'brightness_curves.png', dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), layout='constrained')
    for incidence in (0., 75.):
        for pol, style in (('p', '-'), ('s', '--')):
            d, center = series(rows, incidence, pol, 3, 'cap_dsc_um2')
            _, side = series(rows, incidence, pol, 1, 'cap_dsc_um2')
            _, cp = series(rows, incidence, pol, 3, 'expected_peak_signal_dn')
            _, sp = series(rows, incidence, pol, 1, 'expected_peak_signal_dn')
            axes[0].semilogy(d, center/side, style, marker='o', label='{} deg {}'.format(incidence, pol.upper()))
            axes[1].semilogy(d, cp/sp, style, marker='o', label='{} deg {}'.format(incidence, pol.upper()))
    for ax in axes:
        ax.axhline(1., color='0.5', ls=':')
        ax.set(xlabel='Sphere diameter (nm)', ylabel='CH3 / CH1', title='Collection ratio' if ax is axes[0] else 'Expected camera peak ratio')
        ax.grid(alpha=.2)
        ax.legend(frameon=False)
    fig.suptitle('Central channel is not universally brightest | fixed incident irradiance across angle controls')
    fig.savefig(base/'channel_ratios.png', dpi=180)
    plt.close(fig)
    diameters = (100, 80, 60, 50, 30, 10)
    images = []
    for d in diameters:
        with np.load(base/nist_name/'cases'/('i0_d{:03d}_n64_l12'.format(d))/'raw_fields.npz', allow_pickle=False) as raw:
            images.append([raw['ch{}_p_gray_first_sample_dn'.format(ch)] for ch in (1, 2, 3)])
    low, high = 60., max(float(a.max()) for row in images for a in row)
    fig, axes = plt.subplots(3, len(diameters), figsize=(12, 6), layout='constrained')
    for col, (diameter, row) in enumerate(zip(diameters, images)):
        for ch, array in enumerate(row):
            im = axes[ch, col].imshow(array, origin='lower', cmap='gray', vmin=low, vmax=high, interpolation='nearest')
            axes[ch, col].set_xticks([])
            axes[ch, col].set_yticks([])
            if ch == 0:
                axes[ch, col].set_title('{} nm'.format(diameter))
            if col == 0:
                axes[ch, col].set_ylabel('CH{}'.format(ch+1))
    fig.colorbar(im, ax=list(axes.flat), label='Simulated raw DN, common display window', shrink=.8)
    fig.suptitle('Sphere, 266 nm, normal incidence | same camera settings | no haze / specular leakage / TDI')
    fig.savefig(base/'camera_size_montage.png', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', required=True, type=Path)
    parser.add_argument('--nist-name', default='nist')
    args = parser.parse_args()
    base = args.base
    nist = read(base/args.nist_name/'summary.json')
    if nist['sphere_geometry']['gap_nm'] != 0:
        raise ValueError('This report only accepts spheres touching the wafer')
    convergence, fem = comparisons(base, nist, args.nist_name)
    figures(base, nist, args.nist_name)
    statistics = []
    for incidence in (0., 75.):
        for pol in ('p', 's'):
            for ch in (1, 2, 3):
                d, power = series(nist['rows'], incidence, pol, ch, 'cap_dsc_um2')
                statistics.append(dict(incidence_deg=incidence, polarization=pol, channel=ch,
                    monotonic_increasing=bool(np.all(np.diff(power) > 0)), power_ratio_100_to_10=float(power[-1]/power[0]),
                    adjacent_exponents=local_exponents(d, power).tolist()))
    mesh_steps = mesh_convergence(base)+mesh_convergence(base, ('fem_order2_r022_pml2_m0875_s1',
        'fem_order2_r022_pml2_m075_s1', 'fem_order2_r022_pml2_m0675_s1'))
    result = dict(convergence=convergence, fem_comparison=fem, size_statistics=statistics,
        mesh_convergence=mesh_steps,
        side_symmetry=side_symmetry(nist['rows']),
        max_pupil_vs_cap_relative_change=max(abs(r['pupil_vs_cap_relative_change']) for r in nist['rows']),
        max_flux_relative_error=max(r['flux_relative_error'] for r in nist['rows']))
    from tools.particle_brightness_phase import phase_comparisons
    result['phase_origin_comparisons'] = phase_comparisons(base, args.nist_name)
    (base/'verification.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    from tools.particle_brightness_note import write_note
    write_note(base, nist, result, args.nist_name)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
