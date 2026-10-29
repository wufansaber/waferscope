"""Reproducible 266-nm sphere-size/camera brightness audit, separate from UI code."""
import argparse
from contextlib import contextmanager
import csv
import json
import os
from pathlib import Path
import shutil
import sys

from tools.particle_brightness_metrics import controls, geometry_diagnostic, basis_check, sphere_case


@contextmanager
def native_log(path):
    sys.stdout.flush()
    original = os.dup(1)
    with path.open('w', encoding='utf-8') as log:
        try:
            os.dup2(log.fileno(), 1)
            yield
        finally:
            sys.stdout.flush()
            os.dup2(original, 1)
            os.close(original)


def worker(args):
    control = controls(args.source_run)
    args.output.mkdir(parents=True, exist_ok=True)
    for name in ('request.json', 'imaging_report.json', 'gray_acquisition.png'):
        shutil.copy2(args.source_run/name, args.output/('source_'+name))
    (args.output/'controls.json').write_text(json.dumps(control, indent=2), encoding='utf-8')
    geometry = geometry_diagnostic(control)
    with native_log(args.output/'native_solver.log'):
        checks = dict(basis=basis_check(), geometry=geometry)
        (args.output/'sanity.json').write_text(json.dumps(checks, indent=2), encoding='utf-8')
        rows = []
        for incidence in args.incidences:
            for diameter in range(10, 101, 10):
                folder = args.output/'cases'/('i{:g}_d{:03d}_n64_l12'.format(incidence, diameter))
                rows.extend(sphere_case(control, diameter, incidence, 64, 12, folder))
                (args.output/'progress.json').write_text(json.dumps(dict(
                    completed_cases=len(rows)//6, last_case=str(folder), rows=rows), indent=2), encoding='utf-8')
                print('Completed incidence={} deg, sphere={} nm'.format(incidence, diameter), file=sys.stderr, flush=True)
        for diameter in (() if args.skip_convergence else (10, 100)):
            for samples, lmax in ((96, 12), (64, 16)):
                folder = args.output/'convergence'/('i0_d{:03d}_n{}_l{}'.format(diameter, samples, lmax))
                sphere_case(control, diameter, 0., samples, lmax, folder)
                print('Convergence {}'.format(folder.name), file=sys.stderr, flush=True)
    with (args.output/'brightness.csv').open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = dict(controls=control, sanity=checks, rows=rows,
        sphere_geometry=dict(gap_nm=0., bottom_nm=0., center_height='diameter_nm/2',
                             overrides_archived_gap_nm=control['parameters']['gap_nm']),
        scope='Sphere diameter 10..100 nm, wavelength 266 nm. P/S inputs from the same Jones solution.',
        incidence_control='Incident irradiance is held fixed between incidence-angle controls, not wafer footprint power.',
        limitations=['The source screenshot used a line defect and normal incidence, not a 100-nm sphere.',
                     'Bare-wafer specular illumination and haze are excluded; not a real normal-incidence darkfield SNR prediction.',
                     'BV is an independent scattering reference, but optical propagation/camera reuse the application code.',
                     'Material constants are inherited controlled inputs, not measured UV dispersion.'])
    (args.output/'summary.json').write_text(json.dumps(report, indent=2), encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-run', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--incidences', type=float, nargs='+', default=(0., 75.))
    parser.add_argument('--worker-python', help='Optional separate SCATMECH interpreter; parent provides psutil guard')
    parser.add_argument('--skip-convergence', action='store_true')
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        parser.error('Use an empty output folder to preserve original data')
    if args.worker:
        worker(args)
    else:
        from fullwave_prototype.guarded_run import run_guarded
        result = run_guarded([args.worker_python or sys.executable, '-m', 'tools.particle_brightness_validation', '--worker']+sys.argv[1:],
            args.output.with_suffix('.resources.json'), maximum=6., reserve=3., timeout=3600.,
            cancel_file=args.output.with_suffix('.cancel'))
        if result['status'] != 'completed':
            raise SystemExit(1)


if __name__ == '__main__':
    main()
