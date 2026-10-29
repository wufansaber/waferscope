"""Phase-origin comparison removes one constant phase, never amplitude or position."""
import numpy as np


def constant_phase_relative_l2(actual, reference):
    actual, reference = np.asarray(actual), np.asarray(reference)
    if actual.shape != reference.shape or not np.isfinite(actual).all() or not np.isfinite(reference).all():
        raise ValueError('Matching finite complex fields required')
    norm = np.linalg.norm(reference)
    if norm == 0:
        raise ValueError('Nonzero reference field required')
    phase = np.exp(-1j*np.angle(np.vdot(reference, actual)))
    return float(np.linalg.norm(actual*phase-reference)/norm)


def phase_comparisons(base, nist_name):
    rows = []
    for diameter in (10, 50, 100):
        case = 'i0_d{:03d}_n64_l12'.format(diameter)
        paths = (base/'fem'/('d{}'.format(diameter))/'image_fields.npz',
                 base/'nist'/'cases'/case/'raw_fields.npz',
                 base/nist_name/'cases'/case/'raw_fields.npz')
        if not all(path.is_file() for path in paths) or nist_name == 'nist':
            continue
        with np.load(paths[0], allow_pickle=False) as fem, np.load(paths[1], allow_pickle=False) as old, \
                np.load(paths[2], allow_pickle=False) as corrected:
            if 'ch1_sphere_to_wafer_phase' in old.files:
                continue
            for channel in (1, 2, 3):
                a = fem['ch{}_pupil_with'.format(channel)]
                b = old['ch{}_pupil_modes'.format(channel)][..., 1]
                c = corrected['ch{}_pupil_modes'.format(channel)][..., 1]
                rows.append(dict(diameter_nm=diameter, channel=channel,
                    old_relative_l2=constant_phase_relative_l2(a, b),
                    wafer_reference_relative_l2=constant_phase_relative_l2(a, c),
                    pupil_power_relative_change=float(np.sum(abs(c)**2)/np.sum(abs(b)**2)-1),
                    allowed_alignment='One unobservable constant phase only; no fitted scale or position'))
    return rows
