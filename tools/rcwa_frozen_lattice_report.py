"""Separate lattice angular sampling from periodic-current changes, without refitting."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

import numpy as np

from rcwa_prototype.current_radiation import load_source, lattice_diagnostic


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--lattice', required=True, action='append', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        parser.error('Use an empty output directory to preserve raw data')
    if not args.worker:
        from fullwave_prototype.guarded_run import run_guarded
        result = run_guarded([sys.executable, '-m', 'tools.rcwa_frozen_lattice_report', '--worker']+sys.argv[1:],
            args.output.with_suffix('.resources.json'), maximum=8., reserve=3., timeout=600.)
        if result['status'] != 'completed':
            raise SystemExit(1)
        return
    source, scene, _, provenance = load_source(args.source)
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, path in enumerate(args.lattice):
        grid, lattice_scene, _, grid_provenance = load_source(path)
        for key, value in asdict(scene).items():
            if key not in ('lx_um', 'harmonics_x') and asdict(lattice_scene)[key] != value:
                raise ValueError('Only X lattice and its truncation may change in this diagnostic')
        result = lattice_diagnostic(source, scene, grid, lattice_scene)
        filename = 'lattice_{}.npz'.format(i)
        np.savez_compressed(args.output/filename, **result)
        rows.append(dict(lx_um=lattice_scene.lx_um, raw=filename, lattice=grid_provenance,
                         collected_orders=result['collected_orders'].tolist(),
                         cross_section_um2=result['collected_cross_section_um2'].tolist(),
                         hemisphere_cross_section_um2=float(result['hemisphere_cross_section_um2'])))
    first = np.asarray(rows[0]['cross_section_um2'])
    for row in rows:
        row['relative_change_from_first'] = (np.asarray(row['cross_section_um2'])/first-1).tolist()
    report = dict(source=provenance, scene=asdict(scene), reference_halfwidth=40, lattices=rows,
        status='SAMPLING_ONLY_FROZEN_CURRENT_DIAGNOSTIC',
        limitations=['Both lattices sample the identical selected-cell source; neither is a fresh isolated solve.',
                     'Lattice solid-angle weights are not equal-area cells; coarse cone boundaries can alias.',
                     'No conclusion about convergence of the self-consistent current across cells.'])
    (args.output/'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(rows, indent=2))


if __name__ == '__main__':
    main()
