"""A new large mother requires fresh local FEM; never splice a new draw into old fields."""
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    source = Path(args.source)
    inputs = json.loads((source/'hybrid_input.json').read_text(encoding='utf-8'))
    inputs['mother_range_factor'] = 8
    inputs['coverage'] = .99
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(inputs, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
