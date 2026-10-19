"""SNR from unscaled ADC frame sums, never from an optical photon budget."""
import numpy as np


def measure_roi_series(sample, reference, roi_pixels, rail_values=0):
    sample, reference = np.asarray(sample, float), np.asarray(reference, float)
    if (sample.ndim != 1 or sample.shape != reference.shape or len(sample) < 2
            or not np.isfinite(sample).all() or not np.isfinite(reference).all()):
        raise ValueError('Two or more finite paired ROI frame sums required')
    difference = sample-reference
    mean, noise = float(difference.mean()), float(difference.std(ddof=1))
    valid = rail_values == 0 and noise > 0
    snr = abs(mean)/noise if valid else None
    return dict(frames=len(sample), roi_pixels=int(roi_pixels),
        sample_roi_mean_dn=float(sample.mean()), reference_roi_mean_dn=float(reference.mean()),
        mean_contrast_roi_dn=mean, noise_roi_std_dn=noise, snr=snr,
        nsr=noise/abs(mean) if valid and mean != 0 else None,
        snr_db=float(20*np.log10(snr)) if snr is not None and snr > 0 else None,
        adc_rail_values=int(rail_values), valid=valid,
        definition='abs(mean(sum_ROI(sample_DN - mean_reference_DN))) / '
                   'std(sum_ROI(sample_DN - mean_reference_DN), ddof=1). '
                   'Fixed ROI; raw ADC values, not a fitted template or photon budget.')
