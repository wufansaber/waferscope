"""Local NGSolve field checkpoints for replay without another factorization.

Pickle is used only for locally generated, trusted checkpoints, never external
model uploads. The hash detects accidental changes, not malicious replacement.
"""
from dataclasses import asdict
import hashlib
import json
import pickle
from pathlib import Path


def save(folder, mesh, field, p, geometry, solver):
    import ngsolve as ng
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    data = pickle.dumps((mesh, field), protocol=4)
    (folder/'field.ngs').write_bytes(data)
    metadata = dict(parameters=asdict(p), geometry=asdict(geometry), solver=solver,
        ngsolve_version=ng.__version__, sha256=hashlib.sha256(data).hexdigest(),
        trust='Local workflow checkpoint only; do not load an untrusted pickle.')
    (folder/'checkpoint.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    return metadata


def load_local(folder):
    import ngsolve as ng
    folder = Path(folder)
    metadata = json.loads((folder/'checkpoint.json').read_text(encoding='utf-8'))
    data = (folder/'field.ngs').read_bytes()
    if hashlib.sha256(data).hexdigest() != metadata['sha256']:
        raise ValueError('Checkpoint hash mismatch')
    if ng.__version__ != metadata['ngsolve_version']:
        raise ValueError('Checkpoint requires its original NGSolve version')
    mesh, field = pickle.loads(data)
    return mesh, field, metadata
