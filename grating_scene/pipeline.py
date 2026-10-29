"""Retain native sample/reference records and publish a complete scene atomically."""
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np

from s4_prototype.solver import ROOT, solve
from .collection import collect
from .config import SceneConfig
from .diagnostics import scene_diagnostics
from .imaging import form_images
from .near import air_modal_closure, observation_planes, unpack_planes
from .runtime import default_s4

DEFAULT_EXECUTABLE = default_s4()
LIMITATIONS = [
    'Infinite periodic grating with a particle repeated in each supercell, not an isolated defect.',
    'Plane-wave illumination; power over a reference wafer area is a normalization, not a finite beam.',
    'Native S4 near fields include incident and scattered fields; near maps show |E/E0|^2, not a flux.',
    'Far-field output consists of discrete periodic orders, not a continuous DSC distribution.',
    'Ideal wafer-conjugate Jones relay with paraxial image space, not a full tilted lens/camera model.',
    'Spatial cutoff, particle slicing and isolation must be checked separately; no quantitative camera SNR.',
    'Material indices are explicit user inputs; legacy defaults are not calibrated 266 nm PSL/wafer constants.',
    'This deterministic scene does not include roughness/haze, coating stacks, or detector noise.',
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def is_within(path, folder):
    """Path containment compatible with the application's Python 3.8 environment."""
    try:
        Path(path).resolve().relative_to(Path(folder).resolve())
        return True
    except ValueError:
        return False


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')
    temporary.replace(path)


def _save_arrays(folder, name, arrays):
    if not all(np.isfinite(a).all() for key, a in arrays.items() if key != 'directions'):
        raise ValueError('Nonfinite derived scene data: '+name)
    np.savez_compressed(folder/name, **arrays)


def run_scene(config, output, executable=DEFAULT_EXECUTABLE, timeout=1800):
    config.validate()
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Use an empty scene folder; original calculations are never overwritten')
    output.mkdir(parents=True, exist_ok=True)
    write_json(output/'config.json', config.to_dict())
    planes = observation_planes(config)
    points = np.concatenate([p.reshape(-1, 3) for p in planes.values()]).tolist() if planes else []
    data, native = {}, {}
    for label, particle in (('reference', False), ('sample', config.particle)):
        print('GRATING_STAGE '+label, flush=True)
        data[label], native[label] = solve(config.scene, output/label, particle=particle,
            polarization=config.polarization, near_points=points, executable=executable,
            stripe_rule=True, timeout=timeout)
    return _complete(config, data, native, planes, output)


def _complete(config, data, native, planes, output, reused_from=None):
    print('GRATING_STAGE near_far_collection', flush=True)
    near = unpack_planes(data['reference'], data['sample'], planes)
    air_checks = {key: air_modal_closure(value, config) for key, value in data.items()}
    if any(row.get('sampled') and not row['passed'] for row in air_checks.values()):
        raise RuntimeError('Native near fields do not match reflected modal E/H in exterior air')
    far, collection = collect(data['reference'], data['sample'], config)
    diagnostics = scene_diagnostics(far, native, config)
    if not all(diagnostics['fields'][key]['passive_flux_check_passed'] for key in ('reference', 'sample')):
        raise RuntimeError('Passive scene has negative reflected/transmitted/absorbed flux beyond tolerance')
    images, imaging = form_images(data['reference'], data['sample'], far, config)
    for filename, arrays in (('near_fields.npz', near), ('far_field.npz', far), ('images.npz', images)):
        _save_arrays(output, filename, arrays)
    print('GRATING_STAGE figures', flush=True)
    from .plotting import save_figures
    save_figures(config, near, far, images, collection, output)
    report = dict(format='waferscope-periodic-grating-scene', version=1,
        status='COMPLETED_EXPLORATORY_PERIODIC_SCENE', config=config.to_dict(),
        collection=collection, imaging=imaging, native=native, diagnostics=diagnostics, limitations=LIMITATIONS,
        reused_from=reused_from,
        validation=dict(native_modal_power_checks_passed=True,
            near_fields_sampled=config.near_enabled,
            air_near_far_closure=air_checks,
            image_parseval_checks_passed=True, spatial_cutoff_converged=False,
            isolated_defect_certified=False, real_camera_model_validated=False))
    # The manifest is written last; a cancelled or partial scene has no completion marker.
    report['files_sha256'] = {str(path.relative_to(output)).replace('\\', '/'): digest(path)
        for path in sorted(output.rglob('*')) if path.is_file() and path.suffix in ('.npz', '.lua', '.json', '.txt')}
    write_json(output/'scene_report.json', report)
    print('GRATING_STAGE complete', flush=True)
    return report


def reimage_scene(config, source, output):
    config.validate()
    original, report, _ = load_scene(source)
    frozen = ('scene', 'particle', 'polarization', 'near_enabled', 'near_samples', 'near_xy_z_um', 'near_air_height_um')
    if any(getattr(original, key) != getattr(config, key) for key in frozen):
        raise ValueError('Geometry, illumination polarization or near-plane changes require a new S4 calculation')
    output, source = Path(output), Path(source).resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Use an empty scene folder; original calculations are never overwritten')
    output.mkdir(parents=True, exist_ok=True)
    write_json(output/'config.json', config.to_dict())
    data = {}
    for label in ('reference', 'sample'):
        shutil.copytree(source/label, output/label)
        with np.load(output/label/'modes.npz', allow_pickle=False) as raw:
            data[label] = dict(raw)
    return _complete(config, data, report['native'], observation_planes(config), output,
                     reused_from=str(source))


def load_scene(folder):
    folder = Path(folder).resolve()
    report = json.loads((folder/'scene_report.json').read_text(encoding='utf-8'))
    if (report.get('format') != 'waferscope-periodic-grating-scene' or report.get('version') != 1
            or report.get('status') != 'COMPLETED_EXPLORATORY_PERIODIC_SCENE'):
        raise ValueError('Select a completed periodic grating scene')
    required = {'config.json', 'sample/modes.npz', 'reference/modes.npz',
                'sample/input.lua', 'reference/input.lua', 'near_fields.npz', 'far_field.npz', 'images.npz'}
    if not required.issubset(report['files_sha256']):
        raise ValueError('Scene manifest does not include original native data')
    for relative, expected in report['files_sha256'].items():
        path = (folder/relative).resolve()
        if not is_within(path, folder) or digest(path) != expected:
            raise ValueError('Scene file path/hash check failed: '+relative)
    config = SceneConfig.from_dict(report['config'])
    if SceneConfig.from_dict(json.loads((folder/'config.json').read_text(encoding='utf-8'))) != config:
        raise ValueError('Configuration snapshot does not match the completed report')
    arrays = {}
    for name, filename in (('near', 'near_fields.npz'), ('far', 'far_field.npz'), ('images', 'images.npz')):
        with np.load(folder/filename, allow_pickle=False) as raw:
            arrays[name] = dict(raw)
    return config, report, arrays
