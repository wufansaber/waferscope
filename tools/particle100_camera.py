"""Independent flux-to-electron budgets and Monte Carlo camera readouts."""
from dataclasses import asdict

import numpy as np
from numpy.polynomial.legendre import leggauss

from defect_detection import calculate_detection
from detector_model import integrate_pixels
from fullwave_prototype.direct_imaging import direct_image


def optical_channel(arrays, report, index, name, transmission=1):
    prefix = 'ch{}_'.format(index+1)
    return dict(image_axis_um=arrays[prefix+'image_axis_um'],
        irradiance_w_um2=arrays[prefix+'image_'+name+'_w_um2']*transmission,
        power_w=report['channels'][index][name]['image_power_w']*transmission)


def camera_case(arrays, report, parameters, transmission=1):
    seeds = np.random.SeedSequence(parameters.seed).spawn(3)
    channels = [calculate_detection(optical_channel(arrays, report, i, 'with', transmission),
        optical_channel(arrays, report, i, 'without', transmission),
        report['scattering_parameters']['wavelength_nm'], parameters,
        np.random.default_rng(seeds[i]), sensor_index=i) for i in range(3)]
    return dict(parameters=asdict(parameters), extra_transmission=transmission,
        channels=[{k: v for k, v in c.items() if not isinstance(v, np.ndarray)} for c in channels]), channels


def independent_budget(detector, p, wavelength_nm):
    photon_energy = 6.62607015e-34*299792458/(wavelength_nm*1e-9)
    signal = detector['pixel_power_w']*p.exposure_ms*.001/photon_energy*detector['qe_map']
    haze = detector['haze_signal_e']
    delta = signal-haze
    b = p.exposure_ms*.001*(detector['dark_current_map_e_s']+p.background_e_s)
    variance = (signal+haze/p.reference_frames+(1+1/p.reference_frames)*
        (b+p.read_noise_e**2+p.gain_e_dn**2/12)+(p.haze_mismatch_fraction*haze)**2)
    peak_error = np.max(abs(variance-detector['detection_variance_e2']))/max(np.max(variance), 1e-300)
    template = float(np.sqrt(np.sum(delta**2/variance)))
    return dict(variance_relative_max_error=float(peak_error),
        roi_snr=float(abs(delta.sum())/np.sqrt(variance.sum())), template_snr=template,
        template_relative_error=abs(template/detector['detection_template_snr']-1)
            if detector['detection_template_snr'] else None)


def monte_carlo(detector, p, frames=3000, seed=781):
    """No production camera sampling helpers; per-frame ADC before subtraction."""
    if not detector['detection_roi_valid'] or p.haze_mismatch_fraction != 0:
        raise ValueError('Monte Carlo check needs an unclipped, same-surface reference case')
    rng = np.random.default_rng(seed)
    signal, haze = detector['mean_signal_e'], detector['haze_signal_e']
    background = detector['mean_background_e']
    shape = (frames,)+signal.shape
    def draw(mean):
        charge = rng.poisson(np.broadcast_to(mean, shape))
        readout = np.minimum(charge, p.full_well_e)+rng.normal(0, p.read_noise_e, shape)
        readout += detector['dsnu_offset_e']
        return np.clip(np.rint(readout/p.gain_e_dn+p.offset_dn), 0, 2**p.bit_depth-1)
    observed = draw(signal+background)
    reference = sum(draw(haze+background) for _ in range(p.reference_frames))/p.reference_frames
    difference = (observed-reference)*p.gain_e_dn
    variance, delta = detector['detection_variance_e2'], detector['defect_contrast_e']
    integrated = difference.sum(axis=(1, 2))
    weight = delta/variance
    matched = np.einsum('nij,ij->n', difference, weight)
    roi_snr = abs(integrated.mean())/integrated.std(ddof=1)
    matched_snr = abs(matched.mean())/matched.std(ddof=1)
    result = dict(frames=frames, empirical_roi_snr=float(roi_snr),
        theoretical_roi_snr=detector['detection_roi_snr'],
        empirical_template_snr=float(matched_snr), theoretical_template_snr=detector['detection_template_snr'],
        roi_variance_relative_error=float(abs(integrated.var(ddof=1)/variance.sum()-1)),
        template_snr_relative_error=float(abs(matched_snr/detector['detection_template_snr']-1)))
    result['passed'] = max(result['roi_variance_relative_error'], result['template_snr_relative_error']) < .08
    return result


def direct_pixels(arrays, report, detector, index, parameters):
    """Selected pixel powers: independent vector diffraction plus Gauss pixel integral."""
    prefix, s = 'ch{}_'.format(index+1), report['settings']
    if s['analyzer'] != 'none' or any(s.get('aberrations_nm', (0, 0, 0))):
        raise ValueError('This selected-pixel audit is restricted to the unaberrated, unanalyzed benchmark')
    axis = arrays[prefix+'pupil_u_axis']
    uu, vv = np.meshgrid(axis, axis)
    mask = arrays[prefix+'pupil_mask']
    uv = np.column_stack((uu[mask], vv[mask]))
    weights = np.full(len(uv), (axis[1]-axis[0])**2)
    nx, ny = detector['pixel_power_w'].shape[1], detector['pixel_power_w'].shape[0]
    pixels = [(ny//2-1, nx//2-1), (ny//2-1, nx//2), (ny//2, nx//2-1),
              (ny//2, nx//2), (ny//2+2, nx//2), (ny//2, nx//2+2)]
    half = parameters.pixel_pitch_um*np.sqrt(parameters.fill_factor)/2
    results = {}
    for name in ('with', 'without'):
        integrals = []
        for order in (6, 10):
            nodes, w = leggauss(order)
            points = []
            for row, col in pixels:
                x, y = np.meshgrid(detector['pixel_x_um'][col]+half*nodes,
                                    detector['pixel_y_um'][row]+half*nodes)
                points.append(np.column_stack((x.ravel(), y.ravel())))
            _, _, flux = direct_image(arrays[prefix+'amplitude_'+name][mask], uv, weights,
                np.concatenate(points), report['scattering_parameters']['wavelength_nm']/1000,
                s['tilt_degrees'][index], s['magnification'], s['image_defocus_um'], s['transmission'])
            power = np.sum(flux.reshape(len(pixels), order, order)*np.outer(w, w), axis=(1, 2))
            integrals.append(power*half**2*s['irradiance_w_um2'])
        if name == 'with':
            actual = np.array([detector['pixel_power_w'][r, c] for r, c in pixels])
        else:
            _, _, powers = integrate_pixels(arrays[prefix+'image_axis_um'],
                arrays[prefix+'image_without_w_um2'], parameters)
            actual = np.array([powers[r, c] for r, c in pixels])
        norm = max(np.linalg.norm(integrals[-1]), 1e-300)
        results[name] = dict(pixel_power_relative_l2=float(np.linalg.norm(actual-integrals[-1])/norm),
            oracle_quadrature_relative_l2=float(np.linalg.norm(integrals[0]-integrals[1])/norm))
    return dict(pixels_yx=pixels, components=results,
                scope='Selected six pixels only; same sampled complex pupil, no image alignment or intensity fit.')
