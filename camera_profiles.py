"""Manual-sourced detector presets; no TDI accumulation or calibrated-gain claim."""
from dataclasses import asdict

from detector_model import DetectorParameters


GEMI_NOTE = ('Gemi 8K area approximation, no TDI accumulation. QE >=63.9% at 266 nm; '
    '12.1 e- median read noise and >=16000 e- well from Appendix 3 are assumed area inputs. '
    '16 e-/DN and black level 64 DN are uncalibrated assumptions. '
    'The manual does not guarantee Area-mode image quality or specify whether read noise includes ADC quantization.')


def gemi_area_parameters(gain_e_dn=16., dark_current_e_s=1000.):
    return DetectorParameters(enabled=True, pixel_pitch_um=5., pixels_x=20, pixels_y=20,
        exposure_ms=.01, quantum_efficiency=.639, read_noise_e=12.1,
        dark_current_e_s=dark_current_e_s, full_well_e=16000., gain_e_dn=gain_e_dn,
        bit_depth=10, offset_dn=64., prnu_fraction=.00124, dsnu_e=5.6,
        reference_frames=1, haze_mismatch_fraction=0.)


def describe_profile(name):
    if name not in ('custom', 'gemi_area'):
        raise ValueError('Unknown camera profile')
    return GEMI_NOTE if name == 'gemi_area' else 'Custom detector parameters; calibration is user supplied.'


def profile_values(name):
    if name != 'gemi_area':
        raise ValueError('Choose Gemi area before applying its parameters')
    return asdict(gemi_area_parameters())
