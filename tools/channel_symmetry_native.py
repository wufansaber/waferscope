"""Independently remesh and solve the smooth wafer and mirrored rough wafer."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np

from fullwave_prototype.hybrid_run import load_sources
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.rough_surface import Roughness, build_rough_mesh, height_mesh
from fullwave_prototype.solver import solve
from fullwave_prototype.volume_sources import polarization_sources
from .channel_symmetry import REFLECTION


def extract(mesh, p, enabled):
    _, field, metrics = solve(p, mesh=mesh, particle_enabled=enabled)
    points, weights, volumes = polarization_sources(mesh, field, p, enabled, order=5,
                                                    vectorized_background=True)
    metrics.update(contrast_volumes_um3=volumes, source_quadrature_points=len(points),
                   source_quadrature_order=5)
    return points, weights, metrics


def run(source, output):
    data, metadata = load_sources(source)
    p = Parameters(**metadata['parameters'])
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    report = dict(source=str(Path(source).resolve()), parameters=asdict(p),
                  physical_signature=metadata['physical_signature'],
                  scope='Fresh volume meshes and Maxwell solves at the current preview accuracy; '
                        'not an h/p convergence study. Only the local patch is FEM.')
    print('Stage 1: smooth-wafer Gaussian-beam FEM...', flush=True)
    from fullwave_prototype.solver import build_mesh
    points, weights, metrics = extract(build_mesh(p), p, True)
    np.savez_compressed(output/'flat_sources.npz', source_points=points,
                        source_polarization_weights=weights)
    report['flat'] = metrics
    (output/'native_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    del points, weights

    rough = Roughness(**metadata['roughness'])
    vertices, triangles = height_mesh(rough)
    mirrored_vertices = vertices*REFLECTION
    mirrored_triangles = (triangles*REFLECTION)[:, [0, 2, 1]]
    print('Stage 2: exact mirrored triangular surface, independent mesh...', flush=True)
    # Keep triangle diagonals under reflection; regenerating node heights alone changes them.
    with patch('fullwave_prototype.rough_surface.height_mesh',
               return_value=(mirrored_vertices, mirrored_triangles)):
        mesh = build_rough_mesh(p, rough)
    result = {}
    report['mirror'] = {}
    for enabled, name in ((False, 'without'), (True, 'with')):
        print('Fresh mirrored FEM solve: '+name, flush=True)
        points, weights, metrics = extract(mesh, p, enabled)
        result['source_points_'+name] = points
        result['source_polarization_weights_'+name] = weights
        report['mirror'][name] = metrics
    result['surface_triangles_um'] = mirrored_triangles
    np.savez_compressed(output/'mirror_sources.npz', **result)
    report['mirror']['geometry_height_max_error_um'] = float(np.max(abs(
        mirrored_triangles[:, [0, 2, 1]]*REFLECTION-data['surface_triangles_um'])))
    (output/'native_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('Native symmetry sources saved to '+str(output.resolve()), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    run(args.source, args.output)
