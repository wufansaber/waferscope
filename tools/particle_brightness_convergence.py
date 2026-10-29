"""Extend a completed contact-sphere audit without overwriting any raw case."""
import argparse
import json
from pathlib import Path

from tools.particle_brightness_metrics import sphere_case


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nist', required=True, type=Path)
    args = parser.parse_args()
    report = json.loads((args.nist/'summary.json').read_text(encoding='utf-8'))
    if report['sphere_geometry']['gap_nm'] != 0:
        raise ValueError('Only contact spheres are accepted')
    for incidence in (0., 75.):
        settings = ((64, 12, 32),) if incidence == 0 else ((96, 12, 24), (64, 16, 24), (64, 12, 32))
        for diameter in (10, 100):
            for samples, lmax, cap_order in settings:
                name = 'i{:g}_d{:03d}_n{}_l{}_q{}'.format(incidence, diameter, samples, lmax, cap_order)
                folder = args.nist/'convergence'/name
                if folder.exists():
                    raise ValueError('Refusing to overwrite '+str(folder))
                sphere_case(report['controls'], diameter, incidence, samples, lmax, folder, cap_order)
                print('Completed '+name, flush=True)


if __name__ == '__main__':
    main()
