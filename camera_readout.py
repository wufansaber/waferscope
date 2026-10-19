"""Monochrome CMOS fixed-pattern response and linear ADC, independent of Qt."""
from pathlib import Path

import numpy as np


def sensor_maps(parameters, shape, sensor_index=0):
    if not isinstance(sensor_index, int) or sensor_index < 0:
        raise ValueError('Sensor index must be a nonnegative integer')
    rng = np.random.default_rng(np.random.SeedSequence([parameters.fixed_pattern_seed, sensor_index]))
    sigma = np.sqrt(np.log1p(parameters.prnu_fraction**2))
    response = np.exp(rng.normal(-sigma*sigma/2, sigma, shape))
    qe = np.minimum(parameters.quantum_efficiency*response, 1)
    dark_offset = rng.normal(0, parameters.dsnu_e, shape)
    dark = np.maximum(parameters.dark_current_e_s + rng.normal(0, parameters.dark_current_std_e_s, shape), 0)
    return qe, dark, dark_offset


def digitize(measured_e, nominal_total_e, parameters, fixed_offset_e=0):
    maximum = 2**parameters.bit_depth-1
    analog = measured_e/parameters.gain_e_dn+parameters.offset_dn
    digital = np.clip(np.rint(analog), 0, maximum).astype(np.uint16)
    nominal = (np.minimum(nominal_total_e, parameters.full_well_e)+fixed_offset_e)/parameters.gain_e_dn+parameters.offset_dn
    valid = ((nominal_total_e < parameters.full_well_e) & (nominal > 0) & (nominal < maximum))
    return dict(digital_dn=digital, nominal_dn=nominal, adc_clipped_low=analog < 0,
                adc_clipped_high=analog > maximum, adc_saturated=digital == maximum,
                digital_snr_valid=valid,
                adc_clipped_low_pixels=int(np.count_nonzero(analog < 0)),
                adc_clipped_high_pixels=int(np.count_nonzero(analog > maximum)),
                adc_saturated_pixels=int(np.count_nonzero(digital == maximum)),
                adc_max_dn=maximum, gain_e_dn=parameters.gain_e_dn, offset_dn=parameters.offset_dn)


def export_raw_camera(directory, channels):
    """Lossless 16-bit TIFF containers; no display stretch or peak normalization."""
    from PIL import Image
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for index, detector in enumerate(channels, 1):
        Image.fromarray(detector['digital_dn']).save(directory/'camera_ch{}.tiff'.format(index))
