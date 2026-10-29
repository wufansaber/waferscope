"""Ideal wafer-conjugate Jones relay of periodic complex diffraction orders.

The image space is paraxial. Tilt sets collection and polarization transport;
coordinates are wafer-conjugate, not an uncorrected tilted detector plane.
"""
import numpy as np

from fullwave_prototype.direct_imaging import rotate_between


def mode_coefficients(data, mask, config, tilt):
    indices = np.flatnonzero(mask)
    rays = data['directions'][indices]
    field = data['electric_reflected_z0'][indices]
    angle = np.deg2rad(tilt)
    axis = np.array([np.sin(angle), 0., np.cos(angle)])
    u = np.array([np.cos(angle), 0., -np.sin(angle)])
    axis_array = np.broadcast_to(axis, rays.shape)
    rotated = rotate_between(field, rays, axis_array)
    jones = np.column_stack((rotated@u, rotated[:, 1]))
    jones *= np.sqrt(rays[:, 2]*config.transmission)[:, None]/config.magnification
    if config.analyzer != 'none':
        angle = np.deg2rad({'x': 0, 'y': 90, '45': 45, '135': 135}[config.analyzer])
        vector = np.array([np.cos(angle), np.sin(angle)])
        jones = (jones@vector)[:, None]*vector
    qz = np.sqrt(1-(1-(rays@axis)**2)/config.magnification**2)
    jones *= np.exp(2j*np.pi/config.scene.wavelength_um*config.defocus_um*(qz-1))[:, None]
    return data['orders'][indices], jones


def synthesize(orders, coefficients, samples):
    """Fourier sum with inversion and x=y=0 at the center, no amplitude fit."""
    grid = np.zeros((samples, samples, 2), complex)
    parity = np.where(np.sum(orders, axis=1) % 2, -1., 1.)
    grid[orders[:, 1] % samples, orders[:, 0] % samples] = coefficients*parity[:, None]
    return np.fft.fft2(grid, axes=(0, 1))


def form_images(reference, sample, collection, config):
    n, s, m = config.image_samples, config.scene, config.magnification
    x, y = (np.arange(n)/n-.5)*s.lx_um*m, (np.arange(n)/n-.5)*s.ly_um*m
    pixel_area = s.lx_um*s.ly_um*m*m/(n*n)
    arrays = dict(image_x_um=x, image_y_um=y)
    metrics = []
    for i, (mask, tilt) in enumerate(zip(collection['channel_masks'], config.tilts_deg), 1):
        fields, expected = {}, {}
        for label, data in (('reference', reference), ('sample', sample)):
            orders, coeff = mode_coefficients(data, mask, config, tilt)
            fields[label] = synthesize(orders, coeff, n)
            expected[label] = float(np.sum(abs(coeff)**2)*s.lx_um*s.ly_um*m*m*config.irradiance_w_um2)
            before_analyzer = float(data['reflectance_by_order'][mask].sum()
                * config.incident_cell_power_w*config.transmission)
            if config.analyzer == 'none' and abs(expected[label]-before_analyzer) > 1e-11*max(before_analyzer, 1e-300):
                raise RuntimeError('Polarization transport changed the collected modal power')
        fields['delta'] = fields['sample']-fields['reference']
        density = {name: np.sum(abs(e)**2, axis=-1)*config.irradiance_w_um2 for name, e in fields.items()}
        density['difference'] = density['sample']-density['reference']
        density['interference'] = density['difference']-density['delta']
        row = dict(channel=i, parseval_relative_error={})
        for label in ('sample', 'reference'):
            total = float(density[label].sum()*pixel_area)
            error = abs(total-expected[label])/max(expected[label], 1e-300)
            if error > 1e-10:
                raise RuntimeError('Periodic coherent relay failed its absolute Parseval power check')
            row[label+'_image_power_w'] = total
            row['parseval_relative_error'][label] = error
        for label, data in fields.items():
            arrays['ch{}_electric_{}'.format(i, label)] = data
        for label, data in density.items():
            arrays['ch{}_{}_w_um2'.format(i, label)] = data
        row['difference_image_power_w'] = float(density['difference'].sum()*pixel_area)
        row['delta_field_image_power_w'] = float(density['delta'].sum()*pixel_area)
        metrics.append(row)
    return arrays, dict(channels=metrics, pixel_area_um2=pixel_area,
        model='ideal_wafer_conjugate_periodic_jones_relay', image_space='paraxial',
        coordinate_scope='Uniform wafer-conjugate magnification; not an uncorrected tilted camera plane.',
        phase_scope='Inverted periodic Fourier sum; common Bloch carrier omitted; relative order phase retained.',
        limitation='No camera SNR or isolated-particle certification; no finite illumination or lens prescription.')
