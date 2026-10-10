"""Reassess an existing intensity ensemble without repeating Maxwell solves."""
import argparse
import json
from pathlib import Path

import numpy as np

from fullwave_prototype.haze_statistics import plot_statistics, statistical_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--na', type=float, nargs=3, default=(.5, .5, .5))
    parser.add_argument('--tilts', type=float, nargs=3, default=(45, -45, 0))
    args = parser.parse_args()
    source, output = Path(args.input), Path(args.output)
    report = json.loads((source/'haze_ensemble_report.json').read_text(encoding='utf-8'))
    maps = []
    for seed in report['seeds']:
        with np.load(source/'samples'/('seed_'+str(seed))/'haze_dsc.npz', allow_pickle=False) as data:
            theta, phi = data['theta'], data['phi']
            maps.append(data['dsc'])
    scale = report['metrics']['comparison_area_um2']/report['metrics']['patch_effective_area_um2']
    statistics = statistical_report(maps, theta, phi, scale, args.na, args.tilts)
    statistics['source'] = str(source.resolve())
    output.mkdir(parents=True, exist_ok=True)
    (output/'sample_statistics.json').write_text(json.dumps(statistics, indent=2), encoding='utf-8')
    plot_statistics(output/'haze_sample_convergence.png', statistics)
    print(json.dumps(statistics['observables'], indent=2))


if __name__ == '__main__':
    main()
