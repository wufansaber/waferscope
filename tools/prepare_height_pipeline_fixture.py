"""A known Gaussian-height geometry to test the measured-input imaging path."""
import argparse
import json
from pathlib import Path

import numpy as np

from fullwave_prototype.rough_surface import Roughness, raw_height_grid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    r = Roughness(model='gaussian', grid_points=9)
    np.savez_compressed(path, height_nm=raw_height_grid(r).T*r.rms_nm,
                        dx_um=.025, dy_um=.025)
    path.with_suffix('.json').write_text(json.dumps(dict(
        origin='Synthetic Gaussian-height fixture, not experimental wafer data',
        parent_rms_nm=2, parent_correlation_length_um=.05,
        axes='Rows Y, columns X; no preprocessing intended for geometry equivalence check')),
        encoding='utf-8')


if __name__ == '__main__':
    main()
