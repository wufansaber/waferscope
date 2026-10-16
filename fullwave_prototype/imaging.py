"""Independent ideal objectives fed by complex FEM volume-source radiation."""
from dataclasses import dataclass, asdict, field

import numpy as np

from optical_propagation import pupil_geometry, collimated_pupil, propagate_pupil
from .far_field import radiation_amplitude
from detector_model import DetectorParameters
from defect_detection import calculate_detection
from gray_acquisition import acquire_gray, validate_acquisition
from camera_profiles import describe_profile
from .wavefront import pupil_phase


@dataclass(frozen=True)
class ImagingSettings:
    numerical_apertures: tuple = (.5, .5, .5)
    tilt_degrees: tuple = (45., -45., 0.)
    focal_length_mm: float = 20.
    magnification: float = 20.
    irradiance_w_um2: float = 1e-6
    illumination_power_w: float = None
    illumination_length_um: float = None
    illumination_width_um: float = None
    illumination_incidence_deg: float = None
    illumination_model: str = 'plane'
    illumination_definition: str = 'wafer_gaussian'
    illumination_normalization_area_um2: float = None
    illumination_fit: str = 'manual'
    minimum_field_illumination: float = .5
    transmission: float = 1.
    analyzer: str = 'none'
    image_defocus_um: float = 0.
    pupil_samples: int = 64
    padding: int = 8
    aberrations_nm: tuple = (0., 0., 0.)
    detector: DetectorParameters = field(default_factory=DetectorParameters)
    camera_profile: str = 'custom'
    gray_snr_pairs: int = 0
    gray_snr_seed: int = 20261016

    def validate(self):
        from .field_coverage import FIT_MODES
        if self.illumination_fit not in FIT_MODES or not np.isfinite(self.minimum_field_illumination) or not .05 <= self.minimum_field_illumination <= .95:
            raise ValueError('Invalid camera illumination coverage policy')
        from .beam_inputs import DEFINITIONS
        if self.illumination_definition not in DEFINITIONS:
            raise ValueError('Invalid saved illumination definition')
        if self.illumination_model not in ('plane', 'gaussian'):
            raise ValueError('Unknown illumination power normalization')
        area = self.illumination_normalization_area_um2
        if self.illumination_model == 'gaussian' and (
                area is None or not np.isfinite(area) or area <= 0):
            raise ValueError('Gaussian illumination requires a positive power normalization area')
        if self.illumination_model == 'gaussian' and self.illumination_power_w is None:
            raise ValueError('Gaussian illumination requires total power and footprint dimensions')
        self.detector.validate()
        describe_profile(self.camera_profile)
        validate_acquisition(self.detector, self.gray_snr_pairs, self.gray_snr_seed)
        if len(self.aberrations_nm) != 3 or not np.isfinite(self.aberrations_nm).all():
            raise ValueError('Three finite Zernike coefficients in nm required.')
        if len(self.numerical_apertures) != 3 or len(self.tilt_degrees) != 3:
            raise ValueError('Exactly three NA and tilt values required')
        for value in (*self.numerical_apertures, *self.tilt_degrees, self.focal_length_mm,
                      self.magnification, self.irradiance_w_um2, self.transmission, self.image_defocus_um):
            if not np.isfinite(value):
                raise ValueError('Non-finite imaging setting')
        if any(not 0 <= v < 1 for v in self.numerical_apertures):
            raise ValueError('NA must be in [0,1); zero disables a channel')
        if any(abs(v) > 89.9 for v in self.tilt_degrees):
            raise ValueError('Tilt must be in [-89.9,89.9] degrees')
        if self.magnification < 1 or self.focal_length_mm <= 0 or self.irradiance_w_um2 <= 0:
            raise ValueError('Positive focal length/irradiance and magnification >= 1 required')
        footprint = (self.illumination_power_w, self.illumination_length_um,
                     self.illumination_width_um, self.illumination_incidence_deg)
        if any(value is not None for value in footprint):
            if any(value is None or not np.isfinite(value) for value in footprint):
                raise ValueError('Illumination power, footprint and incidence must be finite and complete')
            power, length, width, incidence = footprint
            if min(power, length, width) <= 0 or not 0 <= incidence < 90:
                raise ValueError('Positive illumination power/footprint and incidence in [0,90) required')
            expected = (power/area if self.illumination_model == 'gaussian' else
                        power/(length*width*np.cos(np.deg2rad(incidence))))
            if not np.isclose(self.irradiance_w_um2, expected, rtol=1e-10, atol=0):
                raise ValueError('Irradiance conflicts with the saved illumination power normalization')
        if not 0 <= self.transmission <= 1:
            raise ValueError('Power transmission must be in [0,1]')
        if self.pupil_samples not in (32, 64, 96, 128) or self.padding not in (2, 4, 8):
            raise ValueError('Use pupil samples 32/64/96/128 and padding 2/4/8')
        if self.analyzer not in ('none', 'x', 'y', '45', '135'):
            raise ValueError('Unknown analyzer')


def axial_flux(e, h):
    return np.real(e[..., 0]*np.conj(h[..., 1])-e[..., 1]*np.conj(h[..., 0]))


def image_from_amplitude(amplitude, local, global_k, rotation, du, wavelength_um, settings):
    # Shared propagation accepts two independent input modes; this solve has one coherent mode.
    modes = np.zeros(amplitude.shape[:-1]+(2, 2), complex)
    modes[..., :, 0] = amplitude
    pupil = collimated_pupil(modes, global_k, local, rotation)
    collected = float(np.sum(abs(pupil)**2)*du**2)
    if settings.analyzer != 'none':
        angle = np.deg2rad({'x': 0, 'y': 90, '45': 45, '135': 135}[settings.analyzer])
        a = np.array((np.cos(angle), np.sin(angle)))
        projection = np.einsum('c,...ci->...i', a, pupil)
        pupil = a[None, None, :, None]*projection[..., None, :]
    pupil *= np.sqrt(settings.transmission)
    if any(settings.aberrations_nm):
        pupil *= pupil_phase(local, du, wavelength_um, settings.aberrations_nm)[..., None, None]
    transmitted = float(np.sum(abs(pupil)**2)*du**2)
    axis, e, h = propagate_pupil(pupil, local, du, wavelength_um, settings)
    e, h = e[..., 0], h[..., 0]
    flux = axial_flux(e, h)
    integral = float(np.sum(flux)*(axis[1]-axis[0])**2)
    error = abs(integral-transmitted)/max(transmitted, 1e-300)
    return dict(axis_um=axis, pupil_xy=pupil[..., 0], e=e, hz0=h, flux=flux,
                collected_dsc_um2=collected, transmitted_dsc_um2=transmitted,
                image_integral_dsc_um2=integral, flux_relative_error=error)


def calculate_channel(data, p, settings, index, cache=None, amplitude_provider=None):
    if amplitude_provider is not None and cache is not None:
        raise ValueError('Custom complex sources cannot use a local-FEM-only pupil cache')
    na, tilt = settings.numerical_apertures[index], settings.tilt_degrees[index]
    local, k, rotation, mask, du = pupil_geometry(na, tilt, settings.pupil_samples)
    arrays = dict(pupil_mask=mask, global_k=k, pupil_u_axis=local[0, :, 0],
                  pupil_axis_mm=local[0, :, 0]*settings.focal_length_mm)
    metrics = dict(na=na, tilt_deg=tilt, enabled=na > 0)
    saved = cache.load(na, tilt, settings.pupil_samples) if cache else None
    metrics['pupil_cache_hit'] = saved is not None
    outputs = {}
    for name in ('with', 'without'):
        amplitude = np.zeros(mask.shape+(2,), complex)
        if saved is not None:
            amplitude = saved[name]
        elif np.any(mask):
            amplitude[mask] = (amplitude_provider(k[mask], name) if amplitude_provider else
                radiation_amplitude(data['source_points_'+name],
                                    data['source_polarization_weights_'+name], k[mask], p))
        result = image_from_amplitude(amplitude, local, k, rotation, du, p.wavelength_nm/1000, settings)
        outputs[name] = result
        arrays['amplitude_'+name] = amplitude
        arrays['pupil_'+name] = result['pupil_xy']
        arrays['image_e_'+name], arrays['image_hz0_'+name] = result['e'], result['hz0']
        arrays['image_'+name+'_w_um2'] = result['flux']*settings.irradiance_w_um2
        metrics[name] = {key: value for key, value in result.items() if not isinstance(value, np.ndarray)}
        metrics[name]['image_power_w'] = result['image_integral_dsc_um2']*settings.irradiance_w_um2
    delta_e = outputs['with']['e']-outputs['without']['e']
    delta_h = outputs['with']['hz0']-outputs['without']['hz0']
    delta_flux = axial_flux(delta_e, delta_h)
    cross = axial_flux(delta_e, outputs['without']['hz0'])+axial_flux(outputs['without']['e'], delta_h)
    arrays['image_axis_um'] = outputs['with']['axis_um']
    arrays['image_change_w_um2'] = delta_flux*settings.irradiance_w_um2
    arrays['image_interference_w_um2'] = cross*settings.irradiance_w_um2
    arrays['image_difference_w_um2'] = arrays['image_with_w_um2']-arrays['image_without_w_um2']
    flux = outputs['with']['flux']
    residual = flux-outputs['without']['flux']-delta_flux-cross
    metrics['decomposition_relative_max_error'] = float(abs(residual).max()/max(abs(flux).max(), 1e-300))
    edge = np.ones(flux.shape, bool)
    border = max(1, len(flux)//10)
    edge[border:-border, border:-border] = False
    metrics['fft_edge_abs_flux_fraction'] = float(abs(flux[edge]).sum()/max(abs(flux).sum(), 1e-300))
    from .illumination import illumination_rays
    metrics['admits_specular'] = bool(na > 0 and any(
        ray[1]@rotation[:, 2] >= np.sqrt(1-na*na) for ray in illumination_rays(p)))
    metrics['specular_test_scope'] = 'Sampled illumination rays only; finite-beam tails are not bounded by this flag.'
    if not all(np.isfinite(a).all() for a in arrays.values()):
        raise RuntimeError('Non-finite imaging result')
    if cache and saved is None:
        cache.save(na, tilt, settings.pupil_samples, {name: arrays['amplitude_'+name] for name in ('with', 'without')})
    return arrays, metrics


def calculate_images(data, p, settings, retain_complex_fields=True, cache_directory=None,
                     amplitude_provider=None):
    p.validate()
    settings.validate()
    if settings.detector.enabled and settings.camera_profile == 'gemi_area' and p.wavelength_nm != 266:
        raise ValueError('Gemi preset QE is specified at 266 nm; select 266 nm or enter a custom calibrated QE')
    from .field_coverage import illumination_coverage
    coverage = illumination_coverage(p, settings.detector, settings.magnification,
        settings.numerical_apertures, settings.tilt_degrees, settings.minimum_field_illumination)
    if settings.illumination_fit == 'camera' and not coverage['all_channels_covered']:
        raise ValueError('Saved illumination does not cover the camera field; calculate new scattering with camera matching')
    if amplitude_provider is not None and cache_directory is not None:
        raise ValueError('Custom complex sources cannot use a local-FEM-only pupil cache')
    arrays, summaries = {}, []
    from .pupil_cache import PupilCache
    cache = PupilCache(cache_directory, data, p) if cache_directory else None
    seeds = np.random.SeedSequence(settings.detector.seed).spawn(3)
    for index in range(3):
        source_label = 'custom complex-source' if amplitude_provider is not None else 'direct FEM'
        print('Computing {} pupil and image for channel {}...'.format(source_label, index+1), flush=True)
        channel, metrics = calculate_channel(data, p, settings, index, cache, amplitude_provider)
        if settings.detector.enabled:
            detector = calculate_detection(dict(image_axis_um=channel['image_axis_um'],
                irradiance_w_um2=channel['image_with_w_um2'],
                power_w=metrics['with']['image_power_w']),
                dict(image_axis_um=channel['image_axis_um'], irradiance_w_um2=channel['image_without_w_um2'],
                     power_w=metrics['without']['image_power_w']), p.wavelength_nm, settings.detector,
                rng=np.random.default_rng(seeds[index]), sensor_index=index)
            detector.update(acquire_gray(detector, settings.detector, settings.gray_snr_pairs,
                                         settings.gray_snr_seed, index))
            detector['camera_profile_note'] = describe_profile(settings.camera_profile)
            metrics['detector'] = {k: v for k, v in detector.items() if not isinstance(v, np.ndarray)}
            channel.update({'detector_'+k: v for k, v in detector.items() if isinstance(v, np.ndarray)})
        if not retain_complex_fields:
            for key in ('image_e_with', 'image_e_without', 'image_hz0_with', 'image_hz0_without'):
                del channel[key]
        arrays.update({'ch{}_{}'.format(index+1, key): value for key, value in channel.items()})
        summaries.append(metrics)
    return arrays, dict(settings=asdict(settings), channels=summaries,
                       field_coverage=coverage,
                       camera_detection_model_version=1 if settings.detector.enabled else 0,
                       tube_focal_length_mm=settings.focal_length_mm*settings.magnification,
                       image_field_units='Flux normalized, not SI E/H; intensity W/um^2.',
                       focus='Object reference plane passes through (0,0,0), perpendicular to each objective axis.',
                       archive_mode='full' if retain_complex_fields else 'compact',
                       scope='Independent aplanatic air objectives with optional Zernike phase and detector. '
                             'Only scattered light; no specular leakage. Detector uses total scattered with-particle flux. '
                             'FEM validation limitations remain; FFT images are periodic.')
