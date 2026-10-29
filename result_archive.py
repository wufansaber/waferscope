"""Self-contained numerical results with safe JSON metadata, never pickle."""
from dataclasses import asdict, is_dataclass
import hashlib
import json
from pathlib import Path
import tempfile
import os

import numpy as np

METADATA_KEY = '_archive_metadata_json'
VERSION = 1


def save_result_archive(path, result, source=None, legacy_arrays=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = dict(legacy_arrays or {})
    aliases = {id(value): key for key, value in arrays.items()}

    def encode(value, parts):
        if isinstance(value, np.ndarray):
            if value.dtype.hasobject:
                raise ValueError('Object arrays are not permitted in raw archives')
            key = aliases.get(id(value), '__'.join(parts))
            if key == METADATA_KEY or (key in arrays and arrays[key] is not value):
                raise ValueError('Archive array name collision: '+key)
            arrays[key] = value
            aliases[id(value)] = key
            return dict(kind='array', key=key)
        if is_dataclass(value) and not isinstance(value, type):
            value = asdict(value)
        if isinstance(value, dict):
            if not all(isinstance(key, str) for key in value):
                raise TypeError('Result dictionary keys must be strings')
            return dict(kind='dict', items={key: encode(item, parts+[key]) for key, item in value.items()})
        if isinstance(value, (list, tuple)):
            return dict(kind='tuple' if isinstance(value, tuple) else 'list',
                        items=[encode(item, parts+[str(i)]) for i, item in enumerate(value)])
        if isinstance(value, np.generic):
            value = value.item()
        if isinstance(value, Path):
            value = str(value)
        if isinstance(value, complex):
            return dict(kind='complex', real=value.real, imag=value.imag)
        if isinstance(value, float) and not np.isfinite(value):
            return dict(kind='nonfinite', value=str(value))
        if value is None or isinstance(value, (str, bool, int, float)):
            return dict(kind='scalar', value=value)
        raise TypeError('Unsupported raw result value: '+type(value).__name__)

    tree = encode(result, [])
    source_tree = encode(source, ['source']) if source is not None else None
    metadata = dict(format='waferscope-raw-result', version=VERSION, result=tree, source=source_tree,
        arrays={key: dict(shape=list(value.shape), dtype=value.dtype.str) for key, value in arrays.items()},
        conventions=dict(angles='rad unless the key ends in _deg',
            raw_data='Stored before plotting; NaNs outside masks and complex phases are preserved.',
            intensity='Native far field: DSC um^2/sr. Objective intensity is a display alias; '
                      'use sampled_dsc_um2_sr/pupil_energy_density and source for raw values.',
            electric='SCATMECH e_x/e_y/e_z are polarization directions, not complex electric fields.',
            coordinates='Use saved x/y/z, global_phi and coordinate metadata; no implicit axis conversion.'))
    arrays[METADATA_KEY] = np.asarray(json.dumps(metadata, allow_nan=False))
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=path.stem+'_', suffix='.npz', delete=False) as stream:
            temporary = Path(stream.name)
            np.savez_compressed(stream, **arrays)
        verify_result_archive(temporary)
        os.replace(str(temporary), str(path))
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            digest.update(chunk)
    manifest = dict(format=metadata['format'], version=VERSION, archive=path.name,
        sha256=digest.hexdigest(), arrays=metadata['arrays'], conventions=metadata['conventions'])
    path.with_suffix('.manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    return manifest


def verify_result_archive(path):
    with np.load(path, allow_pickle=False) as saved:
        metadata = json.loads(str(saved[METADATA_KEY].item()))
        if metadata['format'] != 'waferscope-raw-result' or metadata['version'] != VERSION:
            raise ValueError('Unsupported raw archive format')
        if set(saved.files) != set(metadata['arrays']) | {METADATA_KEY}:
            raise ValueError('Raw archive array inventory differs')
        for key, specification in metadata['arrays'].items():
            value = saved[key]
            if value.dtype.hasobject or list(value.shape) != specification['shape'] or value.dtype.str != specification['dtype']:
                raise ValueError('Raw archive array specification differs: '+key)
    return metadata


def load_result_archive(path):
    """Return plain dictionaries/arrays without importing model classes."""
    metadata = verify_result_archive(path)
    with np.load(path, allow_pickle=False) as saved:
        def decode(node):
            kind = node['kind']
            if kind == 'array':
                return saved[node['key']]
            if kind == 'dict':
                return {key: decode(value) for key, value in node['items'].items()}
            if kind in ('list', 'tuple'):
                values = [decode(value) for value in node['items']]
                return tuple(values) if kind == 'tuple' else values
            if kind == 'complex':
                return complex(node['real'], node['imag'])
            if kind == 'nonfinite':
                return float(node['value'])
            if kind == 'scalar':
                return node['value']
            raise ValueError('Unknown archive node: '+kind)
        return decode(metadata['result']), decode(metadata['source']) if metadata['source'] else None
