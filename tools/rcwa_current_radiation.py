"""Guarded frozen-current continuous-angle diagnostic; not a UI solver."""
import argparse
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--reference', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--level', action='append', type=int, nargs=3, metavar=('RADIAL', 'AZIMUTHAL', 'REF'))
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--timeout', type=float, default=1800.)
    args = parser.parse_args()
    if args.worker:
        from rcwa_prototype.current_radiation import run
        run(args.source, args.reference, args.output, args.level or [(8, 24, 20), (12, 36, 20), (12, 36, 30)])
    else:
        from fullwave_prototype.guarded_run import run_guarded
        result = run_guarded([sys.executable, '-m', 'tools.rcwa_current_radiation', '--worker']+sys.argv[1:],
            args.output.with_suffix('.resources.json'), maximum=8., reserve=3., timeout=args.timeout)
        if result['status'] != 'completed':
            raise SystemExit(1)


if __name__ == '__main__':
    main()
