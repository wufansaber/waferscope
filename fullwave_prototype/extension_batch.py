"""Serial, resource-guarded release checks; failed cases remain in the manifest."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from .guarded_run import run_guarded


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('build/release_upgrade'))
    parser.add_argument('--retry-failed', action='store_true', help='Retry incomplete sensitivity cases at one thread; retain original records')
    args = parser.parse_args()
    entries = json.loads((args.output/'batch.json').read_text(encoding='utf-8')) if args.retry_failed else []
    names = ('beam', 'baseline', 'domain', 'thickness', 'seed8', 'seed9')
    if args.retry_failed:
        complete = {entry['name'] for entry in entries if entry['resources']['status'] == 'completed'}
        names = tuple(name for name in names if name not in complete)
    for name in names:
        folder = name+'_retry1' if args.retry_failed else name
        case = 'beam' if name == 'beam' else 'sensitivity'
        command = [sys.executable, '-m', 'fullwave_prototype.extension_audit', case,
                   '--output', str(args.output/folder)]
        if case == 'sensitivity':
            command += ['--variant', name]
        if args.retry_failed:
            command += ['--threads', '1']
        resources = run_guarded(command, args.output/(folder+'_resources.json'),
                                maximum=12, reserve=3, timeout=900)
        entries.append(dict(name=name, directory=folder, resources=resources))
        (args.output/'batch.json').write_text(json.dumps(entries, indent=2), encoding='utf-8')
    reports = {}
    for entry in entries:
        file = args.output/entry.get('directory', entry['name'])/'audit.json'
        if entry['resources']['status'] == 'completed' and file.is_file():
            reports[entry['name']] = json.loads(file.read_text(encoding='utf-8'))
    powers = {key: np.array([c['with']['image_power_w'] for c in value['channels']])
              for key, value in reports.items() if key != 'beam'}
    summary = dict(completed=list(reports),
                   powers_w={key: value.tolist() for key, value in powers.items()})
    if 'baseline' in powers:
        summary['relative_change_vs_baseline'] = {
            key: (value/powers['baseline']-1).tolist() for key, value in powers.items()}
    if all(key in powers for key in ('baseline', 'seed8', 'seed9')):
        sample = np.array([powers[key] for key in ('baseline', 'seed8', 'seed9')])
        summary['ensemble'] = dict(seeds=[7, 8, 9], mean_power_w=sample.mean(axis=0).tolist(),
            sample_std_w=sample.std(axis=0, ddof=1).tolist(),
            standard_error_w=(sample.std(axis=0, ddof=1)/np.sqrt(3)).tolist(),
            note='Independent intensities averaged, not complex amplitudes; N=3 is not statistical convergence.')
    (args.output/'sensitivity_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2), flush=True)
    raise SystemExit(any(name not in reports for name in ('beam', 'baseline', 'domain', 'thickness', 'seed8', 'seed9')))


if __name__ == '__main__':
    main()
