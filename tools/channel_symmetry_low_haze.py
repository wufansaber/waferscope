"""Make a lower-RMS control without changing illumination or random-surface seed."""
import argparse
import json
from pathlib import Path

from fullwave_prototype.hybrid_run import load_sources


def make_input(source, output, rms_nm):
    if not 0 < rms_nm <= 2:
        raise ValueError('Use a positive diagnostic RMS no greater than the original 2 nm')
    _, metadata = load_sources(source)
    inputs = json.loads((source/'hybrid_input.json').read_text(encoding='utf-8'))
    imaging = json.loads((source/'imaging_report.json').read_text(encoding='utf-8'))
    inputs['parameters'] = metadata['parameters']
    inputs['imaging_settings'] = imaging['settings']
    inputs['roughness']['rms_nm'] = rms_nm
    inputs['assumptions'].append('Lower RMS diagnostic only, not measured wafer roughness. '
        'Same PSD correlation length, seed, grid, taper, beam and optics. Fresh paired local FEM required.')
    path = output/'low_haze_input.json'
    path.write_text(json.dumps(inputs, indent=2), encoding='utf-8')
    print(str(path.resolve()), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('build/field_coverage_20261016/result'))
    parser.add_argument('--output', type=Path, default=Path('build/channel_symmetry_20261016'))
    parser.add_argument('--rms-nm', type=float, default=.2)
    args = parser.parse_args()
    make_input(args.source, args.output, args.rms_nm)
