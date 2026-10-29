"""Inspect SCATMECH exports using only NumPy and the adjacent archive reader."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

try:
    from result_archive import load_result_archive, METADATA_KEY
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from result_archive import load_result_archive, METADATA_KEY


def read_file(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            digest.update(chunk)
    manifest_path = path.with_suffix('.manifest.json')
    verified = None
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        if manifest['sha256'].lower() != digest.hexdigest():
            raise ValueError('NPZ SHA256 differs from the manifest')
        verified = True
    with np.load(path, allow_pickle=False) as saved:
        modern = METADATA_KEY in saved.files
        inventory = {}
        for key in saved.files:
            if key == METADATA_KEY:
                continue
            array = saved[key]
            inventory[key] = dict(shape=list(array.shape), dtype=array.dtype.str,
                                  complex=bool(np.iscomplexobj(array)))
        if not modern:
            result, source = {key: saved[key] for key in saved.files}, None
    if modern:
        result, source = load_result_archive(path)
    summary = dict(file=str(path.resolve()), sha256=digest.hexdigest(),
                   manifest_verified=verified,
                   format='waferscope-raw-result-v1' if modern else 'legacy-flat-npz',
                   plot_kind=result.get('plot_kind', 'native_far_field_or_legacy'),
                   result_keys=list(result), has_source=source is not None, arrays=inventory)
    return result, source, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('file', type=Path)
    parser.add_argument('--inventory', type=Path, help='Write a UTF-8 JSON field inventory')
    args = parser.parse_args()
    _, _, summary = read_file(args.file)
    report = json.dumps(summary, indent=2, ensure_ascii=True)
    if args.inventory:
        args.inventory.parent.mkdir(parents=True, exist_ok=True)
        args.inventory.write_text(report, encoding='utf-8')
    print(report)


if __name__ == '__main__':
    main()
