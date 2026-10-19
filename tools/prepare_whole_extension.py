"""Prepare an immutable fixed-count extension of an existing intensity ensemble."""
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--samples', type=int, default=128)
    args = parser.parse_args()
    report = json.loads((Path(args.source)/'whole_haze_report.json').read_text(encoding='utf-8'))
    inputs = dict(parameters=report['parameters'], roughness=report['roughness'],
        imaging_settings=report['imaging_settings'], samples=args.samples,
        coverage=report['plan']['target_power_coverage'], target_relative_ci=.05)
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(inputs, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
