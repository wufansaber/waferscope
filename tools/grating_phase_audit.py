"""Audit exterior near/far E/H phase closure without rerunning native S4."""
import argparse
from pathlib import Path

import numpy as np

from grating_scene.near import air_modal_closure
from grating_scene.pipeline import digest, load_scene, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scene', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Use a new audit filename')
    config, _, _ = load_scene(args.scene)
    rows = {}
    for key in ('reference', 'sample'):
        with np.load(args.scene/key/'modes.npz', allow_pickle=False) as raw:
            rows[key] = air_modal_closure(dict(raw), config)
    result = dict(scene=str(args.scene.resolve()), scene_report_sha256=digest(args.scene/'scene_report.json'),
                  checks=rows, passed=all(r.get('sampled') and r['passed'] for r in rows.values()))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_json(args.output, result)
    print(result)
    if not result['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
