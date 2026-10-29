"""Copy reusable physical inputs beside images for portable result folders."""
import hashlib
import json
from pathlib import Path
import shutil

LOCAL_FILES = ('report.json', 'paired_fields.npz')
HYBRID_FILES = ('hybrid_fields.npz', 'hybrid_source_report.json', 'hybrid_input.json',
                'shared_surface.npz', 'shared_surface.json',
                'local_shared_height.npz', 'local_shared_height.json',
                'local_scattering/report.json', 'local_scattering/paired_fields.npz')


def digest_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def snapshot_sources(source, destination, hybrid=False):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    inventory = {}
    for name in HYBRID_FILES if hybrid else LOCAL_FILES:
        original, target = source/name, destination/name
        digest = digest_file(original)
        if target != original:
            if target.exists() and digest_file(target) != digest:
                raise ValueError('Refusing to overwrite a different source: '+str(target))
            if not target.exists():
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(original, target)
            if digest_file(target) != digest:
                raise ValueError('Physical source copy failed verification: '+str(target))
        inventory[name] = dict(bytes=target.stat().st_size, sha256=digest)
    destination.mkdir(parents=True, exist_ok=True)
    (destination/'source_manifest.json').write_text(
        json.dumps(dict(format='waferscope-source-snapshot', files=inventory), indent=2), encoding='utf-8')
    return destination
