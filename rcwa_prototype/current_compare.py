"""Compare frozen-cell radiation on identical continuous angular quadratures."""
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import numpy as np

from .current_radiation import matching_scenes
from .geometry import Scene


def compare_arrays(first, second):
    for data in (first, second):
        directions = np.asarray(data['directions'])
        if (directions.shape != data['solid_angle_weights'].shape+(3,)
                or not np.all(np.isfinite(directions))
                or np.any(abs(np.linalg.norm(directions, axis=-1)-1) > 1e-10)):
            raise ValueError('Finite unit directions matching solid-angle weights required')
    for key in ('directions', 'solid_angle_weights'):
        np.testing.assert_allclose(first[key], second[key], atol=1e-14, rtol=1e-13)
    weight = first['solid_angle_weights']
    a, b = first['far_electric_um'], second['far_electric_um']
    if a.shape != b.shape or a.shape != weight.shape+(3,):
        raise ValueError('Matching channel/ray/Cartesian complex-field arrays required')
    if (not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)) or
            not np.all(np.isfinite(weight)) or np.any(weight < 0)):
        raise ValueError('Finite fields and nonnegative solid-angle weights required')
    norm = np.sum(weight*np.sum(abs(a)**2, axis=-1), axis=-1)
    power = np.sum(weight*np.sum(abs(b)**2, axis=-1), axis=-1)
    error = np.sum(weight*np.sum(abs(b-a)**2, axis=-1), axis=-1)
    return [dict(first_cross_section_um2=float(small), second_cross_section_um2=float(large),
                 relative_collected_power_change=float(large/small-1) if small else None,
                 relative_weighted_complex_field_change=float(np.sqrt(delta/small)) if small else None)
            for small, large, delta in zip(norm, power, error)]


def load(folder):
    folder = Path(folder)
    report = json.loads((folder/'report.json').read_text(encoding='utf-8'))
    level = report['levels'][-1]
    raw = folder/level['raw']
    digest = hashlib.sha256(raw.read_bytes()).hexdigest()
    if level.get('raw_sha256') is not None and level['raw_sha256'] != digest:
        raise ValueError('Radiation raw NPZ hash does not match its report')
    with np.load(raw, allow_pickle=False) as archive:
        data = {key: archive[key] for key in archive.files}
    provenance = dict(folder=str(folder.resolve()), raw=level['raw'],
                      sha256=digest, source=report['source'],
                      integrity='recorded_sha256_verified' if 'raw_sha256' in level else 'legacy_no_embedded_radiation_hash')
    return report, data, provenance


def compare_folders(first, second):
    a, fields_a, source_a = load(first)
    b, fields_b, source_b = load(second)
    scene_a, scene_b = Scene(**a['scene']), Scene(**b['scene'])
    if not matching_scenes(scene_a, replace(scene_b, lx_um=scene_a.lx_um, harmonics_x=scene_a.harmonics_x)):
        raise ValueError('Only X cell size and X harmonic count may differ in this diagnostic')
    if (scene_a.harmonics_x is None or scene_b.harmonics_x is None or
            not np.isclose(scene_a.harmonics_x/scene_a.lx_um, scene_b.harmonics_x/scene_b.lx_um,
                           rtol=5e-14, atol=0.)):
        raise ValueError('Matched physical X spatial-frequency cutoff required')
    if a['na'] != b['na'] or a['tilts_deg'] != b['tilts_deg']:
        raise ValueError('Matching collection cones required')
    if a['levels'][-1]['reference_halfwidth'] != b['levels'][-1]['reference_halfwidth']:
        raise ValueError('Matching reciprocal-background truncation required')
    return dict(status='CELL_SENSITIVITY_NOT_ISOLATED_CONVERGENCE',
                first=source_a, second=source_b, channels=compare_arrays(fields_a, fields_b),
                first_lx_um=scene_a.lx_um, second_lx_um=scene_b.lx_um,
                caveat='Same continuous angles and absolute phase reference, no amplitude/phase fit. '
                       'Two periodic cell currents do not establish an isolated-defect limit.')
