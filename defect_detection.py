"""Paired, same-surface camera contrast with haze and EMVA temporal noise."""
import numpy as np

from camera_readout import digitize
from detector_model import calculate_detector


def reference_average(reference, parameters, rng):
    """Average actual ADC frames, retaining per-frame clipping before averaging."""
    mean = reference['digital_dn'].astype(float)
    valid = reference['digital_snr_valid'].copy()
    for _ in range(parameters.reference_frames-1):
        charge = rng.poisson(reference['mean_total_e'])
        measured = np.minimum(charge, parameters.full_well_e).astype(float)
        measured += rng.normal(0, parameters.read_noise_e, measured.shape)
        measured += reference['dsnu_offset_e']
        digital = digitize(measured, reference['mean_total_e'], parameters, reference['dsnu_offset_e'])
        mean += digital['digital_dn']
        valid &= ~(charge >= parameters.full_well_e) & ~digital['adc_clipped_low'] & ~digital['adc_saturated']
    return mean/parameters.reference_frames, valid


def contrast_metrics(sample, reference, parameters, reference_valid):
    """Conditional known-template SNR; fixed sensor patterns are not temporal noise."""
    delta = sample['mean_signal_e']-reference['mean_signal_e']
    n = parameters.reference_frames
    sample_shot = sample['mean_signal_e']
    reference_shot = reference['mean_signal_e']/n
    camera = (sample['mean_background_e']+parameters.read_noise_e**2+parameters.gain_e_dn**2/12
              +(reference['mean_background_e']+parameters.read_noise_e**2+parameters.gain_e_dn**2/12)/n)
    mismatch = (parameters.haze_mismatch_fraction*reference['mean_signal_e'])**2
    variance = sample_shot+reference_shot+camera+mismatch
    snr = np.divide(delta, np.sqrt(variance), out=np.zeros_like(delta), where=variance > 0)
    valid = sample['digital_snr_valid'] & reference_valid
    all_valid = bool(valid.all())
    roi = float(abs(delta.sum())/np.sqrt(variance.sum())) if variance.sum() > 0 else 0.
    template = float(np.sqrt(np.sum(snr**2)))
    return dict(defect_contrast_e=delta, haze_signal_e=reference['mean_signal_e'],
        detection_noise_rms_e=np.sqrt(variance), detection_variance_e2=variance,
        sample_shot_variance_e2=sample_shot, reference_haze_shot_variance_e2=reference_shot,
        camera_variance_e2=camera, haze_mismatch_variance_e2=mismatch,
        detection_snr=snr, detection_valid=valid,
        defect_contrast_roi_e=float(delta.sum()), haze_roi_e=float(reference['mean_signal_e'].sum()),
        detection_noise_roi_e=float(np.sqrt(variance.sum())),
        detection_roi_snr=roi if all_valid else None,
        detection_template_snr=template if all_valid else None,
        valid_pixel_template_snr=float(np.sqrt(np.sum(snr[valid]**2))),
        detection_valid_pixels=int(valid.sum()), detection_invalid_pixels=int((~valid).sum()),
        detection_roi_valid=all_valid, reference_frames=n,
        haze_mismatch_fraction=parameters.haze_mismatch_fraction,
        sample_shot_roi_variance_e2=float(sample_shot.sum()),
        reference_haze_roi_variance_e2=float(reference_shot.sum()),
        camera_roi_variance_e2=float(camera.sum()), haze_mismatch_roi_variance_e2=float(mismatch.sum()),
        detection_status='valid' if all_valid else 'invalid: clipped or unresolved camera pixels',
        detection_scope='Paired same-surface sample minus mean of independent noise reference frames. '
            'Signed contrast includes particle/haze interference and coupling. '
            'Haze photon shot noise in both frames, dark/background, read noise and ADC quantization included. '
            'Fixed DSNU/DCNU offsets cancel; fixed PRNU is retained in pixel response. '
            'Independent pixels and known contrast template; not an empirical detection probability. '
            'Residual haze mismatch RMS is a user assumption, not inferred from one PSD/FEM realization.')


def calculate_detection(sample_channel, haze_channel, wavelength_nm, parameters, rng=None, sensor_index=0):
    if not np.array_equal(sample_channel['image_axis_um'], haze_channel['image_axis_um']):
        raise ValueError('Sample and haze camera axes must match')
    if rng is None:
        rng = np.random.default_rng(parameters.seed)
    sample = calculate_detector(sample_channel, wavelength_nm, parameters, rng, sensor_index)
    reference = calculate_detector(haze_channel, wavelength_nm, parameters, rng, sensor_index)
    average, valid = reference_average(reference, parameters, rng)
    sample.update(contrast_metrics(sample, reference, parameters, valid))
    sample.update(reference_mean_dn=average, reference_digital_dn=reference['digital_dn'],
        reference_mean_e=reference['mean_signal_e'],
        difference_dn=sample['digital_dn'].astype(float)-average,
        reference_power_w=reference['roi_power_w'])
    return sample
