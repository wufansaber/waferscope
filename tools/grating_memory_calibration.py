"""Audit a conservative native-only planning model against completed S4 jobs."""
import argparse
import json
from pathlib import Path

from grating_scene.pipeline import DEFAULT_EXECUTABLE, digest, write_json
from grating_scene.resources import MODEL, native_estimate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, action='append', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Use a new profile path to retain previous measurements')
    executable_hash = digest(DEFAULT_EXECUTABLE)
    rows = []
    for root in args.root:
        for resource_path in root.rglob('*.resources.json'):
            resource = json.loads(resource_path.read_text(encoding='utf-8'))
            if resource.get('status') != 'completed' or 'peak_accounted_gib' not in resource:
                continue
            command = resource.get('command', [])
            if '--output' not in command:
                continue
            folder = Path(command[command.index('--output')+1])
            path = folder/'scene_report.json'
            image_samples = 0
            if path.is_file():
                report = json.loads(path.read_text(encoding='utf-8'))
                if report.get('reused_from'):
                    continue
                info = report['native']['sample']
                image_samples = report['config']['image_samples']
            else:
                path = folder/'report.json'
                if not path.is_file():
                    continue
                info = json.loads(path.read_text(encoding='utf-8'))
            if (info.get('backend') != 'S4' or info.get('executable_sha256') != executable_hash
                    or info.get('factorization') != 'default' or not info.get('stripe_rule')
                    or info.get('one_dimensional') or info['scene'].get('harmonics_x') is None):
                continue
            modes = info['actual_modes']
            layers = info['scene']['slices']+1 if info['particle'] else 1
            prediction = native_estimate(modes, layers, image_samples)['planning_gib']
            peak = resource['peak_accounted_gib']
            rows.append(dict(resource=str(resource_path.resolve()), resource_sha256=digest(resource_path),
                report=str(path.resolve()), report_sha256=digest(path), modes=modes, layers=layers,
                image_samples=image_samples, peak_accounted_gib=peak, planned_gib=prediction,
                measured_headroom_fraction=prediction/peak-1))
    passed = len(rows) >= 3 and all(r['measured_headroom_fraction'] >= .15 for r in rows)
    result = dict(model=MODEL, executable_sha256=executable_hash,
        status='calibration_checks_passed' if passed else 'calibration_checks_failed', observations=rows,
        minimum_measured_headroom_fraction=min((r['measured_headroom_fraction'] for r in rows), default=None),
        guard_limits_unchanged=dict(process_cap_gib=8., system_reserve_gib=3.),
        limitation='Empirical admission only; not an OS limit or a physics validation. Bounded dimension extrapolation.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_json(args.output, result)
    print(json.dumps({key: value for key, value in result.items() if key != 'observations'}, indent=2))
    print('observations:', len(rows))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
