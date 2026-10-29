"""Controlled sphere-brightness verification using existing optics and camera code."""
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path

import numpy as np


def local_exponents(diameters, values):
    d, y = np.asarray(diameters, float), np.asarray(values, float)
    if d.ndim != 1 or d.shape != y.shape or len(d) < 2 or np.any(np.diff(d) <= 0):
        raise ValueError('Matching, increasing diameter and signal series required')
    if not np.isfinite(d).all() or not np.isfinite(y).all() or np.any(d <= 0) or np.any(y <= 0):
        raise ValueError('Finite positive diameter and signal required')
    return np.diff(np.log(y))/np.diff(np.log(d))


def controls(source):
    source = Path(source)
    report = json.loads((source/'imaging_report.json').read_text(encoding='utf-8'))
    p, s = report['scattering_parameters'], report['settings']
    if any(s[key] != 0 for key in ('image_defocus_um',)) or s['analyzer'] != 'none':
        raise ValueError('This benchmark requires focus and no analyzer')
    if s['transmission'] <= 0:
        raise ValueError('Positive optical transmission required')
    return dict(source_folder=str(source.resolve()), parameters=p, settings=s,
                source_sha256={name: hashlib.sha256((source/name).read_bytes()).hexdigest()
                               for name in ('request.json', 'imaging_report.json', 'gray_acquisition.png')})


def geometry_diagnostic(control):
    from fullwave_prototype.parameters import Parameters
    from fullwave_prototype.defects import metadata
    p = Parameters(**control['parameters'])
    sizes = (10., 50., 100.)
    original = [metadata(replace(p, diameter_nm=d)) for d in sizes]
    spheres = [metadata(replace(p, diameter_nm=d, defect_shape='sphere', gap_nm=0.)) for d in sizes]
    return dict(diameters_nm=sizes, original_shape=p.defect_shape, original_geometry=original,
                original_geometry_changes_with_diameter=any(v != original[0] for v in original[1:]),
                sphere_geometry=spheres, sphere_gap_nm=0.,
                sphere_center_heights_nm=[d/2 for d in sizes])


def basis_check():
    from pySCATMECH.local import Local_BRDF_Model
    from pySCATMECH.mueller import Polarization
    model = Local_BRDF_Model('Bobbert_Vlieger_BRDF_Model', dict(
        **{'lambda': .266}, sphere=complex(1.59), radius=.0005, substrate=complex(1.),
        delta=.01, lmax=8, order=-1, Norm_Inc_Approx=0, improve=3))
    incidence = np.deg2rad(35.)
    n = np.array([.3, .4, np.sqrt(.75)])
    theta, phi = np.arccos(n[2]), np.arctan2(-n[0], n[1])
    jones = np.asarray(model.JonesDSC(incidence, theta, phi, coords='psps'))
    columns = np.sum(abs(jones)**2, axis=0)
    mueller = np.array([model.DSC(incidence, theta, phi, coords='psps', inc=Polarization(pol))
                        for pol in ('s', 'p')])
    k0, epsilon = 2*np.pi/.266, 1.59**2
    factor = k0**4*.0005**6*abs((epsilon-1)/(epsilon+2))**2
    rayleigh = np.array([factor*(1-abs(n@e)**2) for e in
                        (np.array([-1., 0., 0.]), np.array([0., np.cos(incidence), np.sin(incidence)]))])
    np.testing.assert_allclose(columns, mueller, rtol=1e-10, atol=0.)
    np.testing.assert_allclose(columns, rayleigh, rtol=5e-4, atol=0.)
    return dict(auxiliary_diameter_nm=1., scope='Free-space Rayleigh limit and s/p column identification',
                jones_column_dsc_um2_sr=columns.tolist(), mueller_dsc_um2_sr=mueller.tolist(),
                rayleigh_dsc_um2_sr=rayleigh.tolist(), relative_rayleigh_change=(columns/rayleigh-1).tolist())


def cap_samples(model, incidence_deg, na, tilt, order=24):
    from fullwave_prototype.aperture_integrals import cap_quadrature
    rays, weights = cap_quadrature(float(na), float(tilt), order)
    theta = np.arccos(rays[:, 2])
    phi = np.arctan2(-rays[:, 0], rays[:, 1])
    jones = np.array([model.JonesDSC(np.deg2rad(incidence_deg), t, f, coords='psps')
                      for t, f in zip(theta, phi)])
    cross = weights@np.sum(abs(jones)**2, axis=1)
    return dict(cap_directions=rays, cap_weights_sr=weights, cap_jones_um=jones), cross


def camera_metrics(channel, mode, wavelength, control, sensor_index):
    from detector_model import DetectorParameters
    from defect_detection import calculate_detection
    from gray_acquisition import acquire_gray
    s = control['settings']
    detector = DetectorParameters(**s['detector'])
    e, h = channel['image_e_modes'][..., mode], channel['image_hz0_modes'][..., mode]
    flux = np.real(e[..., 0]*h[..., 1].conj()-e[..., 1]*h[..., 0].conj())
    axis = channel['image_axis_um']
    power = float(np.sum(flux)*(axis[1]-axis[0])**2*s['irradiance_w_um2'])
    sample = dict(image_axis_um=axis, irradiance_w_um2=flux*s['irradiance_w_um2'], power_w=power)
    reference = dict(image_axis_um=axis, irradiance_w_um2=np.zeros_like(flux), power_w=0.)
    values = calculate_detection(sample, reference, wavelength, detector,
                                 np.random.default_rng(301+sensor_index+mode*3), sensor_index)
    values.update(acquire_gray(values, detector, 100, s['gray_snr_seed']+mode*17, sensor_index))
    expected_dn = values['mean_signal_e']/detector.gain_e_dn
    metrics = dict(image_power_w=power, image_peak_w_um2=float(np.max(sample['irradiance_w_um2'])),
        camera_signal_roi_e=float(values['mean_signal_e'].sum()),
        expected_peak_signal_dn=float(expected_dn.max()), expected_roi_signal_dn=float(expected_dn.sum()),
        expected_peak_raw_dn=float(expected_dn.max()+detector.offset_dn),
        predicted_roi_snr=values['detection_roi_snr'], predicted_template_snr=values['detection_template_snr'],
        acquired_roi_mean_dn=values['gray_contrast_dn'], acquired_roi_noise_dn=values['gray_noise_dn'],
        acquired_roi_snr=values['gray_snr'], camera_valid=values['gray_valid'],
        saturated_pixels=values['saturated_pixels'])
    arrays = {key: value for key, value in values.items() if isinstance(value, np.ndarray)}
    arrays['expected_signal_dn'] = expected_dn
    arrays['image_irradiance_w_um2'] = sample['irradiance_w_um2']
    return metrics, arrays


def sphere_case(control, diameter, incidence, samples, lmax, output, cap_order=24):
    from scatmech_ui.parameters import SimulationParameters, ObjectiveCollectionParameters
    from scatmech_ui.surface_model import build_scattering_model
    from scatmech_ui.imaging import ImagingParameters, _channel
    s, p = control['settings'], control['parameters']
    gap_nm = 0.
    simulation = SimulationParameters(wavelength_nm=266., particle_diameter_nm=float(diameter),
        particle_n=p['particle_n'], particle_k=p['particle_k'], wafer_n=p['substrate_n'],
        wafer_k=p['substrate_k'], incidence_deg=float(incidence), polarization='unpolarized')
    objective = ObjectiveCollectionParameters(focal_length_mm=s['focal_length_mm'],
        channel_numerical_apertures=tuple(s['numerical_apertures']), tilt_degrees=tuple(s['tilt_degrees']))
    imaging = ImagingParameters(magnification=s['magnification'],
        irradiance_mw_um2=s['irradiance_w_um2']*1000, transmission=s['transmission'],
        pupil_samples=samples, padding=4, multipole_order=lmax)
    bundle = build_scattering_model(simulation)
    bundle.model.setParameters(delta=gap_nm/1000, lmax=lmax, order=-1, Norm_Inc_Approx=0, improve=3)
    arrays, rows = {}, []
    for index, (na, tilt) in enumerate(zip(objective.channel_numerical_apertures, objective.tilt_degrees)):
        channel = _channel(bundle.model, simulation, objective, imaging, na, tilt, None)
        caps, cross = cap_samples(bundle.model, incidence, na, tilt, order=cap_order)
        arrays.update({'ch{}_{}'.format(index+1, key): value for key, value in dict(channel, **caps).items()
                       if isinstance(value, np.ndarray)})
        du = float(np.diff(channel['pupil_axis_mm'])[0]/objective.focal_length_mm)
        for mode, polarization in enumerate(('s', 'p')):
            metrics, raw = camera_metrics(channel, mode, 266., control, index)
            pupil_sigma = float(np.sum(abs(channel['pupil_modes'][..., mode])**2)*du**2/imaging.transmission)
            rows.append(dict(metrics, diameter_nm=float(diameter), incidence_deg=float(incidence),
                polarization=polarization, channel=index+1, na=na, tilt_deg=tilt,
                cap_dsc_um2=float(cross[mode]), pupil_dsc_um2=pupil_sigma,
                pupil_vs_cap_relative_change=float(pupil_sigma/cross[mode]-1),
                flux_relative_error=float(abs(metrics['image_power_w']/(pupil_sigma*imaging.transmission*s['irradiance_w_um2'])-1)),
                pupil_samples=samples, lmax=lmax, cap_order=cap_order))
            arrays.update({'ch{}_{}_{}'.format(index+1, polarization, key): value for key, value in raw.items()})
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output/'raw_fields.npz', **arrays)
    report = dict(simulation=asdict(simulation), objective=asdict(objective), imaging=asdict(imaging),
        gap_nm=gap_nm, center_height_nm=diameter/2, sphere_bottom_nm=0.,
        camera_settings=s['detector'], rows=rows,
        raw_sha256=hashlib.sha256((output/'raw_fields.npz').read_bytes()).hexdigest(),
        scope='Exact sphere-on-plane SCATMECH BV with finite multipole/numerical pupil truncation. '
              'Scattered light only, no haze/specular leakage. Shared software optical propagation and detector.',
        field_reference='Top-wafer plane z=0 for propagated fields; native sphere-centered JonesDSC retained unchanged.',
        material_caveat='Fixed legacy material constants from source UI run; not calibrated UV PSL/Si constants.')
    (output/'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return rows
