"""Independent ideal objectives fed by complex FEM volume-source radiation."""
from dataclasses import dataclass, asdict

import numpy as np

from optical_propagation import pupil_geometry, collimated_pupil, propagate_pupil
from .far_field import radiation_amplitude


@dataclass(frozen=True)
class ImagingSettings:
    numerical_apertures: tuple = (.5, .5, .5)
    tilt_degrees: tuple = (45., -45., 0.)
    focal_length_mm: float = 20.
    magnification: float = 20.
    irradiance_w_um2: float = 1e-6
    transmission: float = 1.
    analyzer: str = 'none'
    image_defocus_um: float = 0.
    pupil_samples: int = 64
    padding: int = 8

    def validate(self):
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
    transmitted = float(np.sum(abs(pupil)**2)*du**2)
    axis, e, h = propagate_pupil(pupil, local, du, wavelength_um, settings)
    e, h = e[..., 0], h[..., 0]
    flux = axial_flux(e, h)
    integral = float(np.sum(flux)*(axis[1]-axis[0])**2)
    error = abs(integral-transmitted)/max(transmitted, 1e-300)
    return dict(axis_um=axis, pupil_xy=pupil[..., 0], e=e, hz0=h, flux=flux,
                collected_dsc_um2=collected, transmitted_dsc_um2=transmitted,
                image_integral_dsc_um2=integral, flux_relative_error=error)


def calculate_channel(data, p, settings, index):
    na, tilt = settings.numerical_apertures[index], settings.tilt_degrees[index]
    local, k, rotation, mask, du = pupil_geometry(na, tilt, settings.pupil_samples)
    arrays = dict(pupil_mask=mask, global_k=k, pupil_u_axis=local[0, :, 0],
                  pupil_axis_mm=local[0, :, 0]*settings.focal_length_mm)
    metrics = dict(na=na, tilt_deg=tilt, enabled=na > 0)
    outputs = {}
    for name in ('with', 'without'):
        amplitude = np.zeros(mask.shape+(2,), complex)
        if np.any(mask):
            amplitude[mask] = radiation_amplitude(data['source_points_'+name],
                data['source_polarization_weights_'+name], k[mask], p)
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
    specular = np.array((0., np.sin(np.deg2rad(p.incidence_deg)), np.cos(np.deg2rad(p.incidence_deg))))
    metrics['admits_specular'] = bool(na > 0 and specular@rotation[:, 2] >= np.sqrt(1-na*na))
    if not all(np.isfinite(a).all() for a in arrays.values()):
        raise RuntimeError('Non-finite imaging result')
    return arrays, metrics


def calculate_images(data, p, settings):
    p.validate()
    settings.validate()
    arrays, summaries = {}, []
    for index in range(3):
        print('Computing direct FEM pupil and image for channel {}...'.format(index+1), flush=True)
        channel, metrics = calculate_channel(data, p, settings, index)
        arrays.update({'ch{}_{}'.format(index+1, key): value for key, value in channel.items()})
        summaries.append(metrics)
    return arrays, dict(settings=asdict(settings), channels=summaries,
                       tube_focal_length_mm=settings.focal_length_mm*settings.magnification,
                       image_field_units='Flux normalized, not SI E/H; intensity W/um^2.',
                       focus='Object reference plane passes through (0,0,0), perpendicular to each objective axis.',
                       scope='Independent ideal aplanatic air objectives; no summation of overlapping channels. '
                             'Only scattered light; no specular leakage, aberrations, detector pixels or noise. '
                             'FEM validation limitations remain; FFT images are periodic.')
