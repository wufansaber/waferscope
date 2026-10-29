"""Launch the two frozen end-to-end smoke modes from the packaged directory."""
import argparse
import json
import os
from pathlib import Path

from fullwave_prototype.guarded_run import run_guarded


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--app', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    app, output = args.app.resolve(), args.output.resolve()
    os.environ['QT_QPA_PLATFORM'] = 'offscreen'
    os.chdir(str(app.parent))
    for flag, folder in (('--headless-test', output/'scatmech'), ('--fem-package-test', output)):
        result = run_guarded([str(app), flag, str(folder)],
            output/('scatmech_resources.json' if flag == '--headless-test' else 'fem_resources.json'),
            maximum=6., reserve=3., timeout=300., cancel_file=output/'cancel')
        if result['status'] != 'completed':
            raise RuntimeError('Frozen application smoke failed: '+flag)
    scat = json.loads((output/'scatmech/raw_export_smoke.json').read_text(encoding='utf-8'))
    fem = json.loads((output/'package_smoke.json').read_text(encoding='utf-8'))
    if scat['sphere_gap_um'] != 0 or fem['gap_nm'] != 0:
        raise AssertionError('Frozen application did not retain contact geometry')
    if not scat['coherent_wafer_reference'] or not fem['zero_rhs_skipped_factorization']:
        raise AssertionError('Frozen phase-reference/zero-RHS validation is incomplete')
    print('Frozen SCATMECH and FEM contact smoke checks passed', flush=True)


if __name__ == '__main__':
    main()
