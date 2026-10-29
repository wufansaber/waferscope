"""Optional, bounded empirical admission model for the pinned native S4 build."""
import json
from pathlib import Path

import numpy as np

from .pipeline import digest

MODEL = 'S4-stripes-default-rectangular-v1'


def native_estimate(modes, layers, image_samples=256):
    for name, value, lower in (('modes', modes, 1), ('layers', layers, 1), ('image_samples', image_samples, 0)):
        if isinstance(value, bool) or not isinstance(value, int) or value < lower:
            raise ValueError('Positive integer native dimensions and nonnegative image samples required: '+name)
    one = 16*(2*modes)**2/2**30
    native = .5+one*(4*layers+28)
    image = .5+image_samples**2*16*40/2**30
    return dict(model=MODEL, modes=modes, layers=layers, planning_gib=max(native, image),
        note='Empirical planning with measured headroom, not a hard bound; OS reserve and process cap unchanged.')


def calibrated_estimate(profile_path, executable, modes, layers, image_samples):
    profile = json.loads(Path(profile_path).read_text(encoding='utf-8'))
    if (profile.get('model') != MODEL or profile.get('status') != 'calibration_checks_passed'
            or profile.get('executable_sha256') != digest(executable)):
        raise ValueError('Memory profile does not match this native S4 executable')
    rows = profile['observations']
    if len(rows) < 3:
        raise ValueError('At least three completed native measurements required')
    for row in rows:
        peak = row['peak_accounted_gib']
        predicted = native_estimate(row['modes'], row['layers'], row['image_samples'])['planning_gib']
        if not np.isfinite(peak) or peak <= 0 or predicted < 1.15*peak:
            raise ValueError('Native memory calibration lost its 15% measured headroom')
    if (modes > 1.1*max(r['modes'] for r in rows) or layers > 1.5*max(r['layers'] for r in rows)
            or modes*modes*layers > 1.25*max(r['modes']**2*r['layers'] for r in rows)):
        raise ValueError('Requested dimensions exceed the calibrated memory scope; use the generic planner')
    estimate = native_estimate(modes, layers, image_samples)
    estimate['profile_sha256'] = digest(profile_path)
    return estimate
