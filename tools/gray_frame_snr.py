"""Empirical paired-frame statistics using only unscaled integer camera DN."""
import numpy as np


def measure_paired_frames(sample_dn, reference_dn, roi=None, bit_depth=16):
    sample, reference = np.asarray(sample_dn), np.asarray(reference_dn)
    if (sample.ndim != 3 or sample.shape != reference.shape or sample.shape[0] < 2 or
            min(sample.shape[1:]) == 0):
        raise ValueError('Matching [frame, Y, X] stacks with at least two frames required')
    if not isinstance(bit_depth, int) or isinstance(bit_depth, bool) or not 1 <= bit_depth <= 16:
        raise ValueError('ADC bit depth must be an integer in [1,16]')
    maximum = 2**bit_depth-1
    if any(a.dtype.kind not in 'ui' or a.min() < 0 or a.max() > maximum for a in (sample, reference)):
        raise ValueError('Use unscaled, unsigned/integer raw DN within the ADC range')
    mask = np.ones(sample.shape[1:], bool) if roi is None else np.asarray(roi, bool)
    if mask.shape != sample.shape[1:] or not mask.any():
        raise ValueError('A nonempty, matching fixed ROI is required')
    # Convert before subtraction to preserve negative differences of uint16 images.
    signal = sample[:, mask].astype(np.float64).sum(axis=1)
    background = reference[:, mask].astype(np.float64).sum(axis=1)
    differences = signal-background
    mean = float(differences.mean())
    noise = float(differences.std(ddof=1))
    rails = sum(int(np.count_nonzero((a[:, mask] == 0) | (a[:, mask] == maximum)))
                for a in (sample, reference))
    valid = rails == 0 and noise > 0
    snr = abs(mean)/noise if valid else None
    nsr = noise/abs(mean) if valid and mean != 0 else None
    result = dict(frames=len(sample), roi_pixels=int(mask.sum()),
        sample_roi_mean_dn=float(signal.mean()), reference_roi_mean_dn=float(background.mean()),
        mean_contrast_roi_dn=mean, noise_roi_std_dn=noise, snr=snr, nsr=nsr,
        snr_db=float(20*np.log10(snr)) if snr is not None and snr > 0 else None,
        adc_rail_values=rails, valid=valid,
        definition='abs(mean(sum_ROI(sample_DN-reference_DN))) / '
                   'std(sum_ROI(sample_DN-reference_DN), ddof=1). '
                   'Fixed ROI, independent paired frames, no photon budget or fitted template. '
                   'ADC rails are checked; full-well saturation below the ADC rail needs acquisition metadata.')
    return result, np.column_stack((signal, background, differences))
