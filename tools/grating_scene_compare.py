"""Read-only fixed-cell comparison of completed native periodic scenes."""
import argparse
from pathlib import Path

from grating_scene.pipeline import write_json
from tools.grating_scene_validate import compare_scenes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--refined', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Use a new comparison file')
    result = compare_scenes(args.baseline, args.refined)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_json(args.output, result)
    for row in result['channels']:
        print(row)
    print('Local step passed:', result['local_far_step_passed'])
    print('Isolated defect certified: False')


if __name__ == '__main__':
    main()
