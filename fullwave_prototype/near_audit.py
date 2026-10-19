"""Single particle-present solve with nested near-field sampling and fixed probes."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path

import numpy as np

from .numerical_settings import numerical_parameters
from .rough_surface import Roughness, build_rough_mesh, height_mesh
from .solver import solve, evaluate_fields
from .near_metrics import clearances, peak


def sample_slice(mesh, field, p, y, z):
    yy, zz = np.meshgrid(y, z)
    points = np.column_stack((np.zeros(yy.size), yy.ravel(), zz.ravel()))
    e, _, _ = evaluate_fields(mesh, field, points, p)
    return points, e


def run(scale, output, pml_mesh_scale=1.25, order=3, pml_order=0, pml_inner_fraction=0):
    import ngsolve as ng
    if not np.isfinite(scale) or not .4 <= scale <= 1:
        raise ValueError('Local mesh scale must be between 0.4 and 1')
    p = replace(numerical_parameters('balanced'), gap_nm=10, particle_mesh_um=.015*scale,
                pml_mesh_scale=pml_mesh_scale, order=order, pml_order=pml_order,
                pml_inner_fraction=pml_inner_fraction)
    p.validate()
    roughness = replace(Roughness(), mesh_um=.025*scale)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    ng.SetNumThreads(p.threads)
    mesh = build_rough_mesh(p, roughness)
    _, field, info = solve(p, mesh=mesh)
    return export_diagnostics(mesh, field, p, roughness, info, output, scale)


def export_diagnostics(mesh, field, p, roughness, info, output, scale):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    vertices, triangles = height_mesh(roughness)
    y = np.linspace(-.6*p.domain_radius_um, .6*p.domain_radius_um, 401)
    z = np.linspace(-.299*p.domain_radius_um, .75*p.domain_radius_um, 401)
    points, electric = sample_slice(mesh, field, p, y, z)
    surface, sphere = clearances(points, p, triangles)
    nested = []
    for step in (4, 2, 1):
        indices = np.arange(len(points)).reshape(len(z), len(y))[::step, ::step].ravel()
        nested.append(dict(samples=len(y[::step]), peak=peak(points[indices], electric[indices])))
    # Fixed physical gap probes and nested local sampling, never sample on a boundary.
    zy = np.linspace(-.015, .015, 241)
    zz = np.linspace(vertices[2, 2, 2]+.0002, p.gap_nm/1000-.0002, 241)
    zoom_points, zoom_e = sample_slice(mesh, field, p, zy, zz)
    zs, zb = clearances(zoom_points, p, triangles)
    zoom = []
    for step in (4, 2, 1):
        indices = np.arange(len(zoom_points)).reshape(len(zz), len(zy))[::step, ::step].ravel()
        for distance_nm in (1, 2):
            valid = np.minimum(zs[indices], zb[indices]) >= distance_nm/1000
            zoom.append(dict(samples=len(zy[::step]), clearance_nm=distance_nm,
                             peak=peak(zoom_points[indices], zoom_e[indices], valid)))
    report = dict(parameters=asdict(p), roughness=asdict(roughness), solver=info,
        local_mesh_scale=scale, nested_slice_peaks=nested, gap_nested_peaks=zoom,
        global_sampled_peak=peak(points, electric),
        bulk_sampled_peak=peak(points, electric, (surface > .01) & (sphere > .005)),
        note='Particle-present only; fixed geometry, no smoothing or normalization. YZ slice peaks are not global 3-D maxima. Near-interface maxima may depend on triangulated corners.')
    np.savez_compressed(output/'fields.npz', points=points, electric=electric,
        surface_distance=surface, sphere_distance=sphere, y=y, z=z,
        zoom_points=zoom_points, zoom_e=zoom_e, zoom_surface_distance=zs, zoom_sphere_distance=zb,
        zoom_y=zy, zoom_z=zz, triangles=triangles)
    (output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scale', type=float, default=1)
    parser.add_argument('--output', required=True)
    parser.add_argument('--pml-mesh-scale', type=float, default=1.25)
    parser.add_argument('--order', type=int, choices=(2, 3), default=3)
    parser.add_argument('--pml-order', type=int, choices=(0, 1, 2, 3), default=0)
    parser.add_argument('--pml-inner-fraction', type=float, default=0)
    args = parser.parse_args()
    run(args.scale, args.output, args.pml_mesh_scale, args.order, args.pml_order, args.pml_inner_fraction)


if __name__ == '__main__':
    main()
