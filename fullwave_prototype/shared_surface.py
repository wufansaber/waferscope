"""One synthetic mother surface for whole-beam scattering and the local FEM crop."""
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .experimental_surface import measurement_report
from .whole_surface import random_surface


def surface_hash(axis, height):
    digest = hashlib.sha256()
    for values in (axis, height):
        digest.update(np.asarray(values, dtype='<f8').tobytes())
    return digest.hexdigest()


def make_mother(roughness, plan, range_factor=2):
    if range_factor not in (2, 4, 8):
        raise ValueError('Mother range factor must be 2, 4 or 8; changing it requires new FEM sources')
    n, spacing = range_factor*plan['grid_points'], plan['spacing_um']
    axis = (np.arange(n)-n//2)*spacing
    height = random_surface(roughness, dict(grid_points=n, spacing_um=spacing), roughness.seed)
    return axis, height


def sample_mother(axis, height, x, y):
    """Rows Y, columns X, with an explicit global origin; no half-pixel recentering."""
    xx, yy = np.meshgrid(x, y, indexing='xy')
    interpolate = RegularGridInterpolator((axis, axis), height, bounds_error=True)
    return interpolate(np.stack((yy, xx), axis=-1))


def save_mother(output, axis, height, roughness):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    digest = surface_hash(axis, height)
    np.savez_compressed(output/'shared_surface.npz', axis_um=axis, height_um=height)
    metadata = dict(kind='synthetic_stationary_mother_surface', sha256=digest,
        axis_order='height[Y,X]', origin='Global X=Y=0; stored explicit axis, not image-center inference',
        model_parameters=roughness.__dict__, scope='Synthetic PSD realization, NOT an experimental measurement',
        spacing_um=float(axis[1]-axis[0]), width_um=float(len(axis)*(axis[1]-axis[0])))
    (output/'shared_surface.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    return metadata


def load_mother(source):
    source = Path(source)
    metadata = json.loads((source/'shared_surface.json').read_text(encoding='utf-8'))
    with np.load(source/'shared_surface.npz', allow_pickle=False) as data:
        axis, height = data['axis_um'], data['height_um']
    if surface_hash(axis, height) != metadata['sha256']:
        raise ValueError('Shared mother surface changed; recompute the paired FEM source')
    return axis, height, metadata


def recover_synthesis_parent(source):
    """Recover the actual earlier FFT parent, never repeat or invent exterior heights."""
    from .rough_surface import Roughness
    axis, height, metadata = load_mother(source)
    if metadata.get('kind') != 'synthetic_stationary_mother_surface':
        raise ValueError('Only a stored synthetic PSD mother has a reproducible synthesis parent')
    rough = Roughness(**metadata['model_parameters'])
    n = len(axis)
    if n % 2 or not np.isclose(axis[n//2], 0):
        raise ValueError('Synthesis recovery requires the saved even, origin-centered mother grid')
    step = float(axis[n//2+1]-axis[n//2])
    parent = random_surface(rough, dict(grid_points=n, spacing_um=step), rough.seed,
                            return_synthesis_parent=True)
    crop = parent[n//2:n//2+n, n//2:n//2+n]
    if not np.allclose(crop, height, rtol=1e-12, atol=1e-15):
        raise ValueError('Reconstructed synthesis parent does not reproduce the saved mother surface')
    parent[n//2:n//2+n, n//2:n//2+n] = height
    parent_axis = (np.arange(2*n)-n)*step
    return parent_axis, parent, dict(mother_sha256=metadata['sha256'],
        parent_sha256=surface_hash(parent_axis, parent),
        width_um=float(2*n*step), grid_points=2*n, spacing_um=step,
        central_saved_mother_preserved_exactly=True,
        scope='Earlier finite FFT synthesis parent recovered from the same seed/PSD; '
              'not periodic repetition, not a new draw, not an infinite-wafer extension')


def local_crop(output, axis, height, roughness):
    output = Path(output)
    nodes = np.linspace(-roughness.width_um/2, roughness.width_um/2, roughness.grid_points)
    sampled = sample_mother(axis, height, nodes, nodes)
    path = output/'local_shared_height.npz'
    np.savez_compressed(path, height_nm=sampled*1000,
                        dx_um=nodes[1]-nodes[0], dy_um=nodes[1]-nodes[0])
    (path.with_suffix('.json')).write_text(json.dumps(dict(kind='synthetic_shared_crop',
        mother_sha256=surface_hash(axis, height), detrend='none', coordinate_origin_um=[0, 0],
        note='Synthetic global surface crop, not an independent random realization or measured wafer'),
        indent=2), encoding='utf-8')
    local = replace(roughness, model='measured_height', measurement_path=str(path.resolve()),
                    detrend='none', sampling_grid_levels=[])
    return replace(local, rms_nm=measurement_report(local)['rms_nm'])
