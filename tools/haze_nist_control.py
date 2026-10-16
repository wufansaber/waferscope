"""NIST PSD, finite coherent-beam covariance and saved independent-surface statistics."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path

import numpy as np
from scipy.stats import t

from fullwave_prototype.far_field import angular_grid
from fullwave_prototype.haze_crosscheck import gaussian_reference, _nist_dsc
from fullwave_prototype.hybrid_run import load_sources, optical_settings
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.planar import background_fields
from fullwave_prototype.run import hemisphere_integral
from fullwave_prototype.surface_beam import surface_fields_order
from fullwave_prototype.surface_scattering import radiation_kernel
from optical_propagation import pupil_geometry
from .haze_covariance_reference import expected_tensor, expected_intensity


def nist_pupil_power(directions, geometry, model, p, power_w):
    from pySCATMECH.mueller import Polarization
    incident = np.deg2rad(p.incidence_deg)
    polarization = Polarization(p.polarization)
    intensity = np.array([float((model.MuellerBRDF(incident, np.arccos(q[2]),
        (np.arctan2(q[1], q[0])-np.pi/2) % (2*np.pi), coords='psps')@polarization)[0])
        *q[2]*power_w for q in directions])
    local, _, _, mask, du = geometry
    return float(np.sum(intensity/local[mask, 2])*du**2)


def run(source, output, rms_nm):
    from pySCATMECH.brdf import BRDF_Model
    data, metadata = load_sources(source)
    ensemble_source = json.loads((output/'report.json').read_text(encoding='utf-8'))
    if ensemble_source['physical_signature'] != metadata['physical_signature']:
        raise ValueError('Use the original verified 128-surface source and --rms-nm; '
                         'a different beam or surface source cannot reuse its ensemble')
    p = Parameters(**metadata['parameters'])
    config = json.loads((source/'hybrid_input.json').read_text(encoding='utf-8'))
    settings = optical_settings(config['imaging_settings'], p)
    rough = dict(metadata['stationary_roughness'], rms_nm=rms_nm)
    plane = replace(p, beam_profile='plane')
    theta, phi, q = angular_grid(2)
    frequency2 = (q[..., 0]**2+(q[..., 1]-np.sin(np.deg2rad(p.incidence_deg)))**2)/(p.wavelength_nm/1000)**2
    psd = np.pi*(rms_nm/1000)**2*rough['correlation_length_um']**2*np.exp(
        -(np.pi*rough['correlation_length_um'])**2*frequency2)
    electric = background_fields(np.zeros(3), plane)[0]
    kernel = np.sum(abs(radiation_kernel(q, plane)@electric)**2, axis=-1)
    nist_unit = _nist_dsc(theta, phi, asdict(plane), {None: 'Unit_PSD_Function'}, 1.)
    nist_area = gaussian_reference(theta, phi, asdict(plane), rough, 1.)
    kernel_error = float(np.linalg.norm(kernel-nist_unit)/np.linalg.norm(nist_unit))
    psd_error = float(np.linalg.norm(kernel*psd-nist_area)/np.linalg.norm(nist_area))
    if max(kernel_error, psd_error) > 1e-9:
        raise ValueError('NIST plane-wave normalization or Gaussian PSD convention mismatch')
    axis = data['whole_axis_um']
    step = float(axis[1]-axis[0])
    background = surface_fields_order(axis, axis, p,
        metadata['plan'].get('surface_beam_order', p.beam_order))[0]
    print('Computing finite coherent-beam PSD covariance (no random sampling)...', flush=True)
    tensor = expected_tensor(background*data['whole_height_window'][..., None], step,
                             rms_nm/1000, rough['correlation_length_um'])
    finite = expected_intensity(tensor, step, q, p)*settings.irradiance_w_um2*settings.transmission
    fraction = metadata['plan']['checked_taper_weighted_power_coverage']
    incident_power = settings.illumination_power_w*fraction
    nist = nist_area/np.cos(np.deg2rad(p.incidence_deg))*incident_power*settings.transmission
    model = BRDF_Model('Microroughness_BRDF_Model', {
        'lambda': p.wavelength_nm/1000, 'substrate': complex(p.substrate_n, p.substrate_k),
        'type': 0, 'psd': {None: 'Gaussian_PSD_Function', 'sigma': rms_nm/1000,
                         'length': rough['correlation_length_um']}})
    powers = []
    for i, (na, tilt) in enumerate(zip(settings.numerical_apertures, settings.tilt_degrees)):
        geometry = pupil_geometry(na, tilt, settings.pupil_samples)
        local, directions, _, mask, du = geometry
        directions = directions[mask]
        finite_power = float(np.sum(expected_intensity(tensor, step, directions, p)/local[mask, 2])
                             *du**2*settings.irradiance_w_um2*settings.transmission)
        nist_power = nist_pupil_power(directions, geometry, model, p,
                                     incident_power*settings.transmission)
        powers.append(dict(channel=i+1, finite_coherent_psd_power_w=finite_power,
            nist_plane_psd_power_w=nist_power, finite_vs_nist_relative=finite_power/nist_power-1))
    with np.load(output/'ensemble.npz') as saved:
        samples = saved['powers_w'][:, :, 0]*(rms_nm/metadata['stationary_roughness']['rms_nm'])**2
    for i in (0, 1):
        mean = float(samples[:, i].mean())
        half = float(t.ppf(.975, len(samples)-1)*samples[:, i].std(ddof=1)/np.sqrt(len(samples)))
        expected = powers[i]['finite_coherent_psd_power_w']
        powers[i].update(scaled_128_surface_mean_w=mean, scaled_128_surface_ci_w=[mean-half, mean+half],
            finite_reference_in_95ci=bool(mean-half <= expected <= mean+half),
            ensemble_vs_finite_relative=mean/expected-1)
    finite_total, nist_total = (hemisphere_integral(values, theta) for values in (finite, nist))
    report = dict(rms_nm=rms_nm, correlation_length_um=rough['correlation_length_um'],
        source=str(source.resolve()), parameters=asdict(p),
        incident_power_in_window_w=incident_power, rough_surface_power_coverage=fraction,
        plane_unit_kernel_relative_l2=kernel_error, plane_gaussian_psd_relative_l2=psd_error,
        finite_hemisphere_power_w=finite_total, nist_hemisphere_power_w=nist_total,
        finite_vs_nist_hemisphere_relative=finite_total/nist_total-1, channels=powers,
        scope='Haze only. NIST Microroughness_BRDF_Model is plane-wave first-order Rayleigh-Rice. '
              'Finite Gaussian reference uses full coherent field covariance, not an incoherent sum '
              'of incident rays. 128 SPM surfaces are scaled by RMS^2 exactly within first-order SPM; '
              'they are not 128 fresh coupled FEM solves. Image-plane coherence is not validated '
              'by far-field intensity alone. No empirical intensity normalization or fitted scale.',
        references=['https://pages.nist.gov/SCATMECH/docs/microroughbrdf.htm',
                    'https://pages.nist.gov/SCATMECH/docs/psd.htm'])
    (output/'nist_haze_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    np.savez_compressed(output/'nist_haze_maps.npz', theta=theta, phi=phi,
                        finite_w_sr=finite, nist_plane_w_sr=nist)
    from fullwave_prototype.angular_plotting import angular_disk_edges, format_angular_disk
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
    x, y = angular_disk_edges(theta, phi)
    maximum = max(float(finite.max()), float(nist.max()))
    for ax, values, title in zip(axes, (nist, finite), ('NIST plane-wave PSD', 'Finite coherent-beam PSD covariance')):
        im = ax.pcolormesh(x, y, values, shading='flat', cmap='viridis', vmin=0, vmax=maximum)
        format_angular_disk(ax)
        ax.set_title(title)
    fig.colorbar(im, ax=list(axes), label='W/sr | same absolute scale')
    fig.suptitle('Haze only | RMS {} nm | 266 nm | equal incident power inside roughness window'.format(rms_nm))
    fig.savefig(output/'nist_haze_comparison.png', dpi=170)
    plt.close(fig)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('build/field_coverage_20261016/result'))
    parser.add_argument('--output', type=Path, default=Path('build/channel_symmetry_20261016'))
    parser.add_argument('--rms-nm', type=float, default=.2)
    args = parser.parse_args()
    run(args.source, args.output, args.rms_nm)
