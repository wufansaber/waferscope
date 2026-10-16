"""Bounded-memory single-area-exposure sampling for empirical gray-frame SNR."""
import numpy as np

from camera_readout import digitize
from gray_statistics import measure_roi_series


MAX_VALUES = 200_000_000
MAX_ARCHIVE_VALUES = 6_000_000


def validate_acquisition(p, frames, seed):
    if (isinstance(frames, bool) or not isinstance(frames, int) or frames not in (0,) and not 100 <= frames <= 20000
            or isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed <= 2**32-1):
        raise ValueError('Gray SNR needs 0 (off) or 100..20000 pairs and an unsigned 32-bit seed')
    if p.enabled and frames*p.pixels_x*p.pixels_y*(1+p.reference_frames) > MAX_VALUES:
        raise ValueError('Gray acquisition is too large; reduce ROI, pair count or reference count')


def acquire_gray(detector, p, frames, seed, sensor_index=0):
    validate_acquisition(p, frames, seed)
    if not frames:
        return {}
    if p.haze_mismatch_fraction:
        return dict(gray_status='Not acquired: set haze mismatch to 0 for same-surface temporal SNR',
                    gray_snr=None, gray_valid=False, gray_frames=0)
    shape = detector['mean_signal_e'].shape
    streams = np.random.SeedSequence([seed, sensor_index]).spawn(2)
    rngs = [np.random.default_rng(s) for s in streams]
    sample_sum, ref_sum = np.empty(frames), np.empty(frames)
    archive = frames*np.prod(shape)*(1+p.reference_frames) <= MAX_ARCHIVE_VALUES
    arrays = {}
    if archive:
        arrays['gray_sample_raw_dn'] = np.empty((frames,)+shape, np.uint16)
        arrays['gray_reference_raw_dn'] = np.empty((frames, p.reference_frames)+shape, np.uint16)
    clips, rails = 0, 0
    batch = min(32, max(1, 1_000_000//np.prod(shape)))

    def draw(mean, count, rng):
        nonlocal clips, rails
        total = mean+detector['mean_background_e']
        charge = rng.poisson(np.broadcast_to(total, (count,)+shape))
        clips += int(np.count_nonzero(charge >= p.full_well_e))
        measured = np.minimum(charge, p.full_well_e).astype(float)
        measured += rng.normal(0, p.read_noise_e, measured.shape)
        measured += detector['dsnu_offset_e']
        raw = digitize(measured, total, p, detector['dsnu_offset_e'])['digital_dn']
        rails += int(np.count_nonzero((raw == 0) | (raw == 2**p.bit_depth-1)))
        return raw

    for start in range(0, frames, batch):
        stop = min(frames, start+batch)
        sample = draw(detector['mean_signal_e'], stop-start, rngs[0])
        reference = np.zeros(sample.shape, float)
        for j in range(p.reference_frames):
            raw = draw(detector['haze_signal_e'], stop-start, rngs[1])
            reference += raw
            if archive:
                arrays['gray_reference_raw_dn'][start:stop, j] = raw
            if start == 0 and j == 0:
                arrays['gray_first_reference_dn'] = raw[0].copy()
        reference /= p.reference_frames
        sample_sum[start:stop] = sample.sum(axis=(1, 2), dtype=float)
        ref_sum[start:stop] = reference.sum(axis=(1, 2), dtype=float)
        if archive:
            arrays['gray_sample_raw_dn'][start:stop] = sample
        if start == 0:
            arrays['gray_first_sample_dn'] = sample[0].copy()
            arrays['gray_first_reference_mean_dn'] = reference[0].copy()
    metric = measure_roi_series(sample_sum, ref_sum, np.prod(shape), rails)
    if clips or not detector['detection_roi_valid']:
        metric.update(valid=False, snr=None, snr_db=None)
    arrays['gray_roi_series_dn'] = np.column_stack((sample_sum, ref_sum, sample_sum-ref_sum))
    return dict(arrays, gray_snr=metric['snr'], gray_snr_db=metric['snr_db'], gray_frames=frames,
        gray_contrast_dn=metric['mean_contrast_roi_dn'], gray_noise_dn=metric['noise_roi_std_dn'],
        gray_valid=metric['valid'], gray_status='Valid simulated raw DN' if metric['valid'] else 'Invalid: clipping or unresolved noise',
        gray_full_well_values=clips, gray_adc_rail_values=rails, gray_seed=seed,
        gray_raw_stacks_saved=bool(archive), gray_batch_pairs=int(batch),
        gray_scope='Simulated single area exposures, no TDI multiplication. Same static surface and fixed sensor patterns; '
                   'includes sample/reference haze shot noise, dark current, read noise and actual ADC rounding. '
                   'Temporal ROI SNR, not spatial haze mismatch or hardware measurement. '
                   'Full raw stacks saved only below the archive size cap; ROI DN series and first frames always saved.')
