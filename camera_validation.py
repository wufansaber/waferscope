"""Synthetic EMVA 1288 core checks and a three-channel particle camera example."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path

import numpy as np

from detector_model import DetectorParameters, calculate_detector
from emva1288_noise import pair_statistics, STANDARD_URL


def pixel_integral_check():
    from scipy.special import erf
    p = DetectorParameters(pixels_x=8, pixels_y=6, pixel_pitch_um=1.3, fill_factor=.64,
                           exposure_ms=2, quantum_efficiency=.6)
    axis = np.linspace(-10, 10, 401)
    xx, yy = np.meshgrid(axis, axis)
    sigma, x0, y0, peak = 1.7, .23, -.31, 1e-14
    image = peak*np.exp(-((xx-x0)**2+(yy-y0)**2)/(2*sigma**2))
    channel = dict(image_axis_um=axis, irradiance_w_um2=image, power_w=2*np.pi*sigma**2*peak)
    result = calculate_detector(channel, 266, p)
    half = p.pixel_pitch_um*np.sqrt(p.fill_factor)/2
    def integral(centers, origin):
        return sigma*np.sqrt(np.pi/2)*(erf((centers+half-origin)/(np.sqrt(2)*sigma))
                                      -erf((centers-half-origin)/(np.sqrt(2)*sigma)))
    expected = peak*np.outer(integral(result['pixel_y_um'], y0), integral(result['pixel_x_um'], x0))
    relative = float(np.linalg.norm(result['pixel_power_w']-expected)/np.linalg.norm(expected))
    photons = expected*.002/(6.62607015e-34*299792458/(266e-9))
    electron_error = float(np.linalg.norm(result['mean_signal_e']-.6*photons)/np.linalg.norm(.6*photons))
    return dict(pixel_power_relative_l2=relative, mean_electron_relative_l2=electron_error,
                image_step_um=.05, fill_factor=.64, passed=max(relative, electron_error) < .001,
                reference='Exact displaced Gaussian rectangular-pixel erf integral, not overlap weights.')


def uniform_channel(photons, p):
    axis = np.arange(-10., 11.)
    photon_energy = 6.62607015e-34*299792458/(266e-9)
    irradiance = photons*photon_energy/(p.pixel_pitch_um**2*p.exposure_ms*1e-3)
    return dict(image_axis_um=axis, irradiance_w_um2=np.full((21, 21), irradiance),
                power_w=irradiance*441)


def flat_check():
    p = DetectorParameters(enabled=True, pixels_x=192, pixels_y=192, pixel_pitch_um=.05,
        gain_e_dn=.5, offset_dn=100, quantum_efficiency=.6, read_noise_e=2,
        dark_current_e_s=0, prnu_fraction=.02, dsnu_e=3, fixed_pattern_seed=4)
    levels = np.linspace(0, 25000, 50)
    pairs = []
    for index, photons in enumerate(levels):
        channel = uniform_channel(photons, p)
        a = calculate_detector(channel, 266, replace(p, seed=1000+index*2))
        b = calculate_detector(channel, 266, replace(p, seed=1001+index*2))
        pairs.append(pair_statistics(a['digital_dn'], b['digital_dn']))
    means = np.array([row['mean_dn'] for row in pairs])
    variances = np.array([row['temporal_variance_dn2'] for row in pairs])
    signal_dn = means-means[0]
    # Weight by the approximate variance uncertainty, not by the image brightness.
    design = np.column_stack((signal_dn, np.ones(len(signal_dn))))
    weights = 1/np.maximum(variances, 1)
    slope, intercept = np.linalg.lstsq(design*weights[:, None], variances*weights, rcond=None)[0]
    responsivity = np.polyfit(levels, signal_dn, 1)[0]
    qe = responsivity/slope
    dark_noise = np.sqrt(max(variances[0]-1/12, 0))/slope
    stacks = {}
    for name, photons in (('dark', 0), ('bright', 12500)):
        stack = np.stack([calculate_detector(uniform_channel(photons, p), 266,
            replace(p, seed=5000+i))['digital_dn'].astype(float) for i in range(64)])
        temporal = float(stack.var(axis=0, ddof=1).mean())
        average = stack.mean(axis=0)
        spatial = float(average.var()-temporal/len(stack))
        stacks[name] = dict(mean_dn=float(average.mean()), spatial_variance_dn2=spatial)
    dsnu = np.sqrt(max(stacks['dark']['spatial_variance_dn2'], 0))/slope
    prnu = np.sqrt(max(stacks['bright']['spatial_variance_dn2']-stacks['dark']['spatial_variance_dn2'], 0))
    prnu /= stacks['bright']['mean_dn']-stacks['dark']['mean_dn']
    result = dict(parameters=asdict(p), exposure_levels=50, stack_frames=64,
        recovered_gain_dn_e=float(slope), gain_relative_error=float(abs(slope*p.gain_e_dn-1)),
        recovered_qe=float(qe), qe_relative_error=float(abs(qe/p.quantum_efficiency-1)),
        recovered_temporal_dark_noise_e=float(dark_noise),
        temporal_dark_noise_relative_error=float(abs(dark_noise/p.read_noise_e-1)),
        recovered_dsnu_e=float(dsnu), dsnu_relative_error=float(abs(dsnu/p.dsnu_e-1)),
        recovered_prnu_fraction=float(prnu), prnu_relative_error=float(abs(prnu/p.prnu_fraction-1)),
        fitted_dark_variance_dn2=float(intercept), scope='Synthetic uniform input, white pixel maps; '
            'EMVA pair/stack equations. No uneven illumination, spatial high-pass compliance '
            'measurement, row/column correlation, temperature sweep or camera certification.')
    result['passed'] = max(result[key] for key in ('gain_relative_error', 'qe_relative_error',
        'temporal_dark_noise_relative_error', 'dsnu_relative_error', 'prnu_relative_error')) < .03
    return result, dict(photons=levels, signal_dn=signal_dn, variance_dn2=variances,
                         fitted_variance_dn2=slope*signal_dn+intercept)


def particle_example(output):
    from scatmech_ui.imaging import ImagingParameters, calculate_imaging, export_imaging
    from scatmech_ui.parameters import SimulationParameters, ObjectiveCollectionParameters
    from camera_validation_plot import example_figures
    simulation = replace(SimulationParameters(), wavelength_nm=266, particle_diameter_nm=10,
                         roughness_enabled=False, illumination_length_um=2, illumination_width_um=2)
    irradiance = 1/(2*2*np.cos(np.deg2rad(simulation.incidence_deg)))
    detector = DetectorParameters(enabled=True, pixels_x=16, pixels_y=16, exposure_ms=.01,
        gain_e_dn=.5, offset_dn=100, quantum_efficiency=.6, prnu_fraction=.01, dsnu_e=.5,
        seed=17, fixed_pattern_seed=3)
    imaging = ImagingParameters(pupil_samples=128, padding=4,
                                irradiance_mw_um2=irradiance*1000, detector=detector)
    result = calculate_imaging(simulation, ObjectiveCollectionParameters(), imaging)
    optical, camera = example_figures(result)
    export_imaging(result, output, optical, camera)
    import matplotlib.pyplot as plt
    plt.close(optical)
    plt.close(camera)
    return dict(scattering=asdict(simulation), camera=asdict(detector),
        illumination='Local uniform plane wave, reference irradiance from 1 W / (2x2 um wafer footprint*cos(incidence)); '
                     'not a solved finite 2x2 um beam. Material n,k are input values, not measured dispersion.',
        channels=[dict(channel=i+1, roi_power_w=c['detector']['roi_power_w'],
            signal_e=c['detector']['signal_e'], total_roi_snr=c['detector']['emva_total_roi_snr'],
            full_well_saturated=c['detector']['saturated_pixels'],
            adc_saturated=c['detector']['adc_saturated_pixels']) for i, c in enumerate(result['channels'])])


def saved_fem_camera(source, output):
    """Re-evaluate electronics from archived irradiance without changing scattering."""
    source, output = Path(source), Path(output)
    report = json.loads((source/'imaging_report.json').read_text(encoding='utf-8'))
    p = DetectorParameters(**report['settings']['detector'])
    seeds = np.random.SeedSequence(p.seed).spawn(3)
    arrays, summaries = {}, []
    with np.load(source/'image_fields.npz') as saved:
        for i in range(3):
            prefix = 'ch{}_'.format(i+1)
            channel = dict(image_axis_um=saved[prefix+'image_axis_um'],
                irradiance_w_um2=saved[prefix+'image_with_w_um2'],
                power_w=report['channels'][i]['with']['image_power_w'])
            detector = calculate_detector(channel, report['scattering_parameters']['wavelength_nm'], p,
                                          np.random.default_rng(seeds[i]), sensor_index=i)
            arrays.update({prefix+'detector_'+k: v for k, v in detector.items() if isinstance(v, np.ndarray)})
            summaries.append({k: v for k, v in detector.items() if not isinstance(v, np.ndarray)})
    output.mkdir(parents=True, exist_ok=True)
    from fullwave_prototype.detector_plot import plot_detector
    plot_detector(output, arrays)
    np.savez_compressed(output/'camera_fields.npz', **arrays)
    payload = dict(source_directory=str(source.resolve()), detector_parameters=asdict(p), channels=summaries,
                   scope='Saved FEM optical intensity; fresh camera electronics. No new FEM solve.')
    (output/'camera_report.json').write_text(json.dumps(payload, indent=2), encoding='utf-8')
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='build/camera_validation_20261012')
    parser.add_argument('--saved-image', help='Optional FEM image directory for electronics-only validation')
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    flat, arrays = flat_check()
    from camera_validation_plot import plot_flat
    plot_flat(output, arrays, flat)
    example = particle_example(output/'particle_10nm')
    pixel = pixel_integral_check()
    payload = dict(standard_url=STANDARD_URL, flat=flat, pixels=pixel, particle_example=example)
    if args.saved_image:
        payload['fem_camera_example'] = saved_fem_camera(args.saved_image, output/'fem_camera')
    (output/'report.json').write_text(json.dumps(payload, indent=2), encoding='utf-8')
    np.savez_compressed(output/'flat_checks.npz', **arrays)
    print(json.dumps(payload, indent=2), flush=True)
    if not flat['passed'] or not pixel['passed']:
        raise RuntimeError('Camera synthetic validation failed; inspect report.json')


if __name__ == '__main__':
    main()
