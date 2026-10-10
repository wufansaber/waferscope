"""Run a native rough-mesh regression in an isolated resource-bounded process."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

from fullwave_prototype.parameters import Parameters
from fullwave_prototype.rough_surface import Roughness, build_rough_mesh, height_mesh


def audit(request, output, particle=False, cad_scale=1000.):
    import ngsolve as ng
    p, r = Parameters(**request['parameters']), Roughness(**request['roughness'])
    ng.SetNumThreads(p.threads)
    started = time.perf_counter()
    print('Building audited rough mesh...', flush=True)
    mesh = build_rough_mesh(p, r, particle_enabled=particle, cad_scale=cad_scale)
    # Geometric volumes must exclude the artificial complex PML coordinate stretch.
    mesh.UnSetPML('pml_.*')
    vertices, triangles = height_mesh(r)
    volumes = {name: float(ng.Integrate(1, mesh, definedon=mesh.Materials(name)))
               for name in set(mesh.GetMaterials())}
    signed_volume = r.width_um**2*triangles[..., 2].mean()
    mesh_signed = volumes.get('bump', 0)-volumes.get('pit', 0)
    expected_total = 4*np.pi*(p.domain_radius_um+p.pml_thickness_um)**3/3
    report = dict(elements=mesh.ne, volumes_um3=volumes,
                  signed_height_volume_um3=float(signed_volume),
                  signed_volume_error_um3=float(mesh_signed-signed_volume),
                  total_volume_relative_error=float(abs(sum(volumes.values())/expected_total-1)),
                  surface_sha256=hashlib.sha256(vertices.tobytes()).hexdigest(),
                  elapsed_seconds=time.perf_counter()-started)
    failures = []
    if abs(mesh_signed-signed_volume) > max(r.width_um**2*r.rms_nm/1000*1e-5, 1e-12):
        failures.append('Mesh does not preserve the signed roughness volume')
    if report['total_volume_relative_error'] > 1e-3:
        failures.append('Mesh regions do not fill the spherical domain')
    report.update(status='failed' if failures else 'passed', failures=failures,
                  cad_scale=cad_scale, volume_measure='Physical geometric volume; PML disabled for this audit')
    (output/'mesh_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)
    if failures:
        raise RuntimeError('; '.join(failures))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--particle', action='store_true')
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--cad-scale', type=float, default=1000.)
    args = parser.parse_args()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.worker:
        audit(json.loads(Path(args.input).read_text(encoding='utf-8')), output, args.particle, args.cad_scale)
        return
    from fullwave_prototype.guarded_run import run_guarded
    command = [sys.executable, '-m', 'tools.rough_mesh_audit', '--worker',
               '--input', args.input, '--output', str(output), '--cad-scale', str(args.cad_scale)]
    if args.particle:
        command.append('--particle')
    result = run_guarded(command, output/'resources.json', maximum=8, timeout=120)
    raise SystemExit(0 if result['status'] == 'completed' else 1)


if __name__ == '__main__':
    main()
