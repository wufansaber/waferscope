"""Validate a fixed-configuration intensity ensemble before adding independent seeds."""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import numpy as np


def load_resume(source, p, roughness, settings, samples):
    source = Path(source).resolve()
    report = json.loads((source/'whole_haze_report.json').read_text(encoding='utf-8'))
    if report.get('status') != 'complete' or report.get('model') != 'whole_beam_first_order_surface':
        raise ValueError('Only a completed whole-beam SPM ensemble can be extended')
    physical, saved_p = asdict(p), dict(report['parameters'])
    from .parameters import Parameters
    saved_p = asdict(Parameters(**saved_p))
    physical.pop('beam_order')
    saved_p.pop('beam_order')
    if physical != saved_p or asdict(roughness) != report['roughness']:
        raise ValueError('Ensemble extension requires identical scattering inputs')
    requested, saved = asdict(settings), dict(report['imaging_settings'])
    saved.setdefault('illumination_definition', 'wafer_gaussian')
    for item in (requested, saved):
        item.pop('irradiance_w_um2')
        item.pop('illumination_normalization_area_um2')
        item.pop('detector')
    if json.loads(json.dumps(requested)) != saved:
        raise ValueError('Ensemble extension requires identical optics and power')
    count = report['metrics']['samples']
    if not isinstance(count, int) or not count < samples <= 128:
        raise ValueError('Choose a fixed final count above the existing count, at most 128')
    if report['seeds'] != list(range(roughness.seed, roughness.seed+count)):
        raise ValueError('Saved seed sequence is not the requested independent sequence')
    archive = source/'whole_haze_fields.npz'
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    if report.get('fields_sha256') and digest != report['fields_sha256']:
        raise ValueError('Saved ensemble fields changed')
    with np.load(archive, allow_pickle=False) as saved:
        arrays = {key: saved[key] for key in saved.files}
    if any(not np.isfinite(values).all() for values in arrays.values()):
        raise ValueError('Saved ensemble contains nonfinite fields')
    observations = np.asarray(report['statistics']['sample_observables_um2'])
    if observations.shape != (count, 4) or np.any(observations < 0) or not np.isfinite(observations).all():
        raise ValueError('Invalid saved sample observables')
    plan = report['plan']
    axis = (np.arange(plan['grid_points'])-plan['grid_points']//2)*plan['spacing_um']
    if not np.allclose(axis, arrays['x_um'], rtol=1e-12, atol=1e-12):
        raise ValueError('Saved surface grid differs from the ensemble plan')
    return report, arrays, dict(source=str(source), previous_samples=count, fixed_final_samples=samples,
        source_fields_sha256=digest, selection='Fixed final count chosen before extension, no repeated-peeking stop',
        old_source_integrity='Verified content hash' if report.get('fields_sha256') else
                             'Legacy source without prior checksum; captured current checksum and checked configuration')
