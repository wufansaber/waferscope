"""Read-only physical consistency audit of a completed periodic scene."""
import argparse

from grating_scene.diagnostics import scene_diagnostics
from grating_scene.pipeline import digest, load_scene, write_json


def main():
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Use a new audit filename; previous records are never overwritten')
    config, report, arrays = load_scene(args.source)
    audit = scene_diagnostics(arrays['far'], report['native'], config)
    audit.update(source=str(args.source.resolve()),
        source_report_sha256=digest(args.source/'scene_report.json'))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_json(args.output, audit)


if __name__ == '__main__':
    main()
