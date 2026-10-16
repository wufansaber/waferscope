"""Whole-illumination coherent first-order haze and independent image planes."""
import argparse
from dataclasses import asdict, replace
import json
import hashlib
from pathlib import Path

import numpy as np

from optical_propagation import pupil_geometry
from detector_model import DetectorParameters
from .far_field import angular_grid
from .finite_illumination import power_area
from .haze_statistics import report_from_observables
from .imaging import ImagingSettings, image_from_amplitude
from .parameters import Parameters
from .rough_surface import Roughness
from .surface_scattering import SurfaceRadiation, direct_amplitude
from .whole_surface import plan_surface, random_surface, window
from .surface_observables import checked_integrals


def run(inputs, output, resume_source=None):
    p = Parameters(**inputs['parameters'])
    p.validate()
    roughness = Roughness(**inputs['roughness'])
    values = dict(inputs['imaging_settings'])
    detector = DetectorParameters(**values.pop('detector', {}))
    if detector.enabled:
        raise ValueError('Whole-beam haze exports optical flux, not detector noise; disable the detector')
    settings = ImagingSettings(**values, detector=detector)
    settings.validate()
    from .field_coverage import illumination_coverage
    field_coverage = illumination_coverage(p, detector, settings.magnification,
        settings.numerical_apertures, settings.tilt_degrees, settings.minimum_field_illumination)
    if settings.illumination_fit == 'camera' and not field_coverage['all_channels_covered']:
        raise ValueError('Beam does not cover the camera field; refit illumination and recompute')
    if (settings.illumination_model != 'gaussian' or settings.illumination_definition != p.beam_definition or
            not np.allclose((settings.illumination_length_um, settings.illumination_width_um,
                             settings.illumination_incidence_deg),
                            (p.beam_length_um, p.beam_width_um, p.incidence_deg), rtol=1e-10, atol=0)):
        raise ValueError('Illumination settings must match the coherent scattering beam')
    samples = inputs.get('samples', 16)
    if not isinstance(samples, int) or not 1 <= samples <= 128:
        raise ValueError('Use 1..128 independent rough surfaces')
    if roughness.seed+samples-1 > 1000000:
        raise ValueError('Last surface seed exceeds 1000000')
    resume = None
    if resume_source:
        from .whole_haze_resume import load_resume
        previous_report, previous_arrays, resume = load_resume(resume_source, p, roughness, settings, samples)
        p, plan = Parameters(**previous_report['parameters']), previous_report['plan']
    else:
        p, plan = plan_surface(p, roughness, inputs.get('coverage', .99),
            inputs.get('variance_fraction', .9), inputs.get('taper_width_um', .025))
    settings = replace(settings, illumination_normalization_area_um2=power_area(p),
                       irradiance_w_um2=settings.illumination_power_w/power_area(p))
    settings.validate()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    n, step, width = plan['grid_points'], plan['spacing_um'], plan['width_um']
    from .field_coverage import rough_surface_coverage
    plan['camera_coverage'] = rough_surface_coverage(settings, plan)
    axis = (np.arange(n)-n//2)*step
    edge = window(axis, axis, width, plan['taper_width_um'])
    from .surface_beam import surface_fields_order
    background, ei, hi = surface_fields_order(axis, axis, p, plan.get('surface_beam_order', p.beam_order))
    flux = -np.real(np.cross(ei, hi.conj()))[..., 2]
    coverage = float(np.sum(flux*edge**2)*step**2/power_area(p))
    if not plan['target_power_coverage'] <= coverage <= 1.01:
        raise ValueError('Fine-grid taper-weighted power coverage failed: {:.6g}'.format(coverage))
    plan['checked_taper_weighted_power_coverage'] = coverage
    plan['checked_geometric_power_coverage'] = float(np.sum(flux)*step**2/power_area(p))
    theta, phi, directions = angular_grid(2.5)
    map_sum, image_sums, channels = np.zeros(directions.shape[:-1]), {}, []
    image_powers = []
    realized = []
    observable_series, integral_checks = [], []
    first_height = None
    transform_checks = []
    start = 0
    if resume:
        start = resume['previous_samples']
        if not np.array_equal(theta, previous_arrays['theta']) or not np.array_equal(phi, previous_arrays['phi']):
            raise ValueError('Saved ensemble angular grid mismatch')
        map_sum = previous_arrays['mean_dsc_um2_sr']*start
        first_height = previous_arrays['first_height_nm']/1000
        channels = previous_report['channels']
        realized = [previous_report['metrics']['stationary_sample_rms_mean_nm']]*start
        image_powers = previous_arrays['sample_image_power_w'].tolist()
        observable_series = previous_report['statistics']['sample_observables_um2'][:]
        integral_checks = previous_report['direct_integration_checks'][:]
        transform_checks = previous_report['sample_selected_fft_amplitude_errors'][:]
        for channel in range(3):
            image_sums[channel] = previous_arrays['ch{}_mean_image_w_um2'.format(channel+1)]*start
            image_sums['axis'+str(channel)] = previous_arrays['ch{}_image_axis_um'.format(channel+1)]
    for index in range(start, samples):
        print('Whole-beam first-order haze surface {}/{} | {} nodes/axis'.format(index+1, samples, n), flush=True)
        height = random_surface(roughness, plan, roughness.seed+index)
        realized.append(float(np.sqrt(np.mean(height**2))*1000))
        if first_height is None:
            first_height = height*edge
        weighted = height[..., None]*edge[..., None]*background
        radiation = SurfaceRadiation(axis, axis, weighted,
                                     p, padding=inputs.get('fft_padding', 4))
        selected = np.array([[.1, .2, np.sqrt(.95)], [-.3, .5, np.sqrt(.66)],
                             [.7, -.2, np.sqrt(.47)], [.1, -.8, np.sqrt(.35)]])
        direct = direct_amplitude(axis, axis, weighted, selected, p)
        fft = radiation.amplitude(selected)
        transform_check = float(np.linalg.norm(fft-direct)/max(np.linalg.norm(direct), 1e-300))
        transform_checks.append(transform_check)
        if transform_check > .01:
            raise ValueError('Surface FFT interpolation differs from direct quadrature by >1%; increase FFT padding')
        amplitude = radiation.amplitude(directions)
        map_sum += np.sum(abs(amplitude)**2, axis=-1)
        observable, check = checked_integrals(radiation, settings.numerical_apertures, settings.tilt_degrees)
        observable_series.append(observable)
        integral_checks.append(check)
        powers = []
        for channel in range(3):
            na, tilt = settings.numerical_apertures[channel], settings.tilt_degrees[channel]
            local, k, rotation, mask, du = pupil_geometry(na, tilt, settings.pupil_samples)
            a = np.zeros(mask.shape+(2,), complex)
            if np.any(mask):
                a[mask] = radiation.amplitude(k[mask])
            result = image_from_amplitude(a, local, k, rotation, du, p.wavelength_nm/1000, settings)
            image = result['flux']*settings.irradiance_w_um2
            image_sums[channel] = image_sums.get(channel, np.zeros_like(image))+image
            image_sums['axis'+str(channel)] = result['axis_um']
            powers.append(result['image_integral_dsc_um2']*settings.irradiance_w_um2)
            if index == 0:
                channels.append(dict(na=na, tilt_deg=tilt, enabled=na > 0,
                    image_flux_relative_error=result['flux_relative_error'],
                    object_fft_field_width_um=float((result['axis_um'][-1]-result['axis_um'][0])/settings.magnification)))
        image_powers.append(powers)
        del radiation
    statistics = report_from_observables(observable_series, settings.numerical_apertures,
        settings.tilt_degrees, inputs.get('target_relative_ci', .05))
    if resume:
        resume['extension_only_statistics'] = report_from_observables(observable_series[start:],
            settings.numerical_apertures, settings.tilt_degrees, inputs.get('target_relative_ci', .05))
    statistics['aperture_integration'] = 'Direct complex surface spectrum on solid-angle quadrature; '
    statistics['aperture_integration'] += 'hemisphere Gauss/uniform azimuth; cone-aligned clipped apertures'
    power_means = {name: item['mean_um2']*settings.irradiance_w_um2
                  for name, item in statistics['observables'].items()}
    warnings = ['First-order smooth-interface approximation, NOT a full-wafer FEM solve. '
                'No multiple roughness scattering, particle coupling or specular reflection.',
                'Finite surface synthesis and FFT images need bandwidth, padding and size convergence.',
                'Only the taper-weighted covered power is modeled; no missing-area extrapolation.']
    if samples < 16:
        warnings.append('Few independent samples: no claim of ensemble convergence.')
    if not statistics['statistical_target_met']:
        warnings.append('The requested statistical confidence half-width has not been met for all active channels.')
    if not all(check['target_met'] for check in integral_checks):
        warnings.append('A direct angular integration failed the 1% adjacent-order target at order 128.')
    if roughness.model == 'abc' and roughness.psd_exponent <= 4:
        warnings.append('ABC C<=4 has no finite infinite-band RMS slope; grid cutoff remains physical uncertainty.')
    for channel in channels:
        if channel['enabled'] and channel['object_fft_field_width_um'] < width:
            warnings.append('An image FFT object field is narrower than the rough surface; increase pupil samples '
                            'before quantitative image-shape use. Padding alone cannot enlarge the field.')
    report = dict(status='complete', model='whole_beam_first_order_surface', source_directory=str(output.resolve()),
        field_coverage=field_coverage,
        parameters=asdict(p), roughness=asdict(roughness), plan=plan,
        imaging_settings=asdict(settings), statistics=statistics,
        direct_integration_checks=integral_checks,
        selected_fft_vs_direct_amplitude_relative_l2=max(transform_checks),
        sample_selected_fft_amplitude_errors=transform_checks,
        seeds=list(range(roughness.seed, roughness.seed+samples)), channels=channels,
        metrics=dict(samples=samples, area_extrapolation_applied=False,
            mean_power_w=power_means, mean_image_power_w=np.mean(image_powers, axis=0).tolist(),
            stationary_sample_rms_mean_nm=float(np.mean(realized))), warnings=warnings,
        definition='Each realization is coherent over the complete illuminated rough surface; '
                   'independent realizations are averaged in intensity, not amplitude.')
    if resume:
        report['extension'] = resume
    arrays = dict(theta=theta, phi=phi, mean_dsc_um2_sr=map_sum/samples,
        x_um=axis, y_um=axis, first_height_nm=first_height*1000,
        incident_flux_w_um2=flux*settings.irradiance_w_um2, height_window=edge,
        sample_image_power_w=np.array(image_powers))
    for channel in range(3):
        arrays['ch{}_mean_image_w_um2'.format(channel+1)] = image_sums[channel]/samples
        arrays['ch{}_image_axis_um'.format(channel+1)] = image_sums['axis'+str(channel)]
    np.savez_compressed(output/'whole_haze_fields.npz', **arrays)
    report['fields_sha256'] = hashlib.sha256((output/'whole_haze_fields.npz').read_bytes()).hexdigest()
    from .whole_haze_plotting import plot_results
    plot_results(output, arrays, report)
    (output/'whole_haze_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report['metrics'], indent=2), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--resume', help='Completed ensemble; extend with a fixed final sample count')
    args = parser.parse_args()
    run(json.loads(Path(args.input).read_text(encoding='utf-8')), args.output, args.resume)


if __name__ == '__main__':
    main()
