"""Same-contact-scene S4/grcwa complex-field check; not a bare-grating oracle."""
import argparse
import json
from pathlib import Path
import shutil
import sys

import numpy as np

from grating_scene.collection import collect
from grating_scene.diagnostics import scene_diagnostics
from grating_scene.pipeline import digest, load_scene, write_json
from rcwa_prototype.geometry import matrix_memory_estimate
from tools.rcwa_crosscheck_report import summarize, write_outputs


def coincident_xy(s4_arrays, grcwa_data):
    points = s4_arrays['xy_points_um']
    if not np.isclose(float(grcwa_data['near_z_um']), points[0, 0, 2], rtol=0, atol=1e-13):
        raise ValueError('Near planes differ; no field fitting/interpolation allowed')
    a, b = [], []
    for j, y in enumerate(points[:, 0, 1]):
        jy = np.flatnonzero(abs(grcwa_data['near_y_um']-y) < 1e-13)
        for i, x in enumerate(points[0, :, 0]):
            ix = np.flatnonzero(abs(grcwa_data['near_x_um']-x) < 1e-13)
            if len(ix) == len(jy) == 1:
                a.append(s4_arrays['xy_sample_electric'][j, i])
                b.append(grcwa_data['near_electric'][jy[0], ix[0]])
    if not a:
        raise ValueError('Near grids have no coincident samples')
    a, b = np.asarray(a), np.asarray(b)
    return dict(points=len(a), relative_complex_e_l2=float(np.linalg.norm(b-a)/max(np.linalg.norm(a), 1e-300)),
        scope='Total E on exactly coincident XY points; no interpolation or fitted amplitude/phase.')


def run(source, output):
    from rcwa_prototype.solver import solve
    config, native, arrays = load_scene(source)
    if not config.particle or config.scene.kind != 'grating' or config.polarization != 'p':
        raise ValueError('This check requires a particle plus grating with p polarization')
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Use an empty crosscheck folder')
    output.mkdir(parents=True, exist_ok=True)
    write_json(output/'config.json', config.to_dict())
    data, info, paths = {}, {}, {}
    near = bool(arrays['near']) and config.near_xy_z_um == .02
    for label, particle in (('reference', False), ('sample', True)):
        print('HYBRID_CROSSCHECK '+label, flush=True)
        target = output/('s4_'+label)
        shutil.copytree(Path(source)/label, target)
        paths['s4_'+label] = str(target.resolve())
        data[label], info[label] = solve(config.scene, particle=particle, near=near, particle_method='analytic')
        target = output/('grcwa_'+label)
        target.mkdir()
        np.savez_compressed(target/'modes.npz', **data[label])
        info[label].update(particle=particle, raw_npz_sha256=digest(target/'modes.npz'))
        write_json(target/'report.json', info[label])
        paths['grcwa_'+label] = str(target.resolve())
    result = summarize(paths)
    far, powers = collect(data['reference'], data['sample'], config)
    result.update(s4_scene_report_sha256=digest(Path(source)/'scene_report.json'),
        contact_gap_um=config.scene.gap_um, grcwa_collection=powers,
        grcwa_diagnostics=scene_diagnostics(far, info, config),
        s4_diagnostics=scene_diagnostics(arrays['far'], native['native'], config))
    if near:
        result['sample_near_xy'] = coincident_xy(arrays['near'], data['sample'])
        reference_near = dict(arrays['near'], xy_sample_electric=arrays['near']['xy_reference_electric'])
        result['reference_near_xy'] = coincident_xy(reference_near, data['reference'])
    result['quantitative_hybrid_certified'] = False
    result['limitations'].append('New contact sphere series; implementation agreement does not certify cutoff or isolation.')
    write_outputs(output, result)
    write_json(output/'crosscheck.json', result)
    print(json.dumps(result['fields'], indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--timeout', type=float, default=3600.)
    args = parser.parse_args()
    config, _, _ = load_scene(args.source)
    if args.worker:
        run(args.source, args.output)
        return
    import psutil
    from fullwave_prototype.guarded_run import run_guarded
    estimate = matrix_memory_estimate(config.scene.requested_modes, config.scene.slices+1)
    available = psutil.virtual_memory().available/2**30
    blocked = estimate['planning_gib'] > min(8., available-3.)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_json(args.output.with_suffix('.admission.json'), dict(planning=estimate, available_gib=available,
        system_reserve_gib=3., maximum_gib=8., config=config.to_dict(),
        status='blocked_memory_preflight' if blocked else 'admitted'))
    if blocked:
        parser.error('Memory preflight blocked this matched mixed-scene check')
    result = run_guarded([sys.executable, '-m', 'tools.grating_hybrid_crosscheck', '--worker']+sys.argv[1:],
        args.output.with_suffix('.resources.json'), maximum=8., reserve=3., timeout=args.timeout)
    if result['status'] != 'completed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
