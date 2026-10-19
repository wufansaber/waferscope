"""Separate h/p check for a symmetry defect found in the smooth-wafer preview."""
import argparse
from dataclasses import asdict, replace
import gc
import json
from pathlib import Path

import numpy as np

from detector_model import DetectorParameters
from fullwave_prototype.hybrid_run import optical_settings
from fullwave_prototype.parameters import Parameters
from optical_propagation import pupil_geometry
from .channel_symmetry import channel_image, relative_error, reflect_centered_grid
from .gemi_area_profile import profile


def solve_refinements(output, fine_order2=False):
    from fullwave_prototype.solver import build_mesh
    from .channel_symmetry_native import extract
    native = json.loads((output/'native_report.json').read_text(encoding='utf-8'))
    baseline = Parameters(**native['parameters'])
    report = dict(source=native['source'], physical_signature=native['physical_signature'], cases={})
    previous_path = output/'flat_refinement_native.json'
    if previous_path.exists():
        previous = json.loads(previous_path.read_text(encoding='utf-8'))
        if previous['physical_signature'] != native['physical_signature']:
            raise ValueError('Refinement cache does not match the current physical source')
        report = previous
    cases = [('order2', 2, 1., False), ('order1_h075', 1, .75, False)]
    if fine_order2:
        cases.append(('order2_h075_condensed', 2, .75, True))
    for name, order, scale, condensed in cases:
        p = replace(baseline, order=order, air_mesh_um=baseline.air_mesh_um*scale,
                    substrate_mesh_um=baseline.substrate_mesh_um*scale,
                    particle_mesh_um=baseline.particle_mesh_um*scale,
                    static_condensation=condensed)
        if (report['cases'].get(name, {}).get('parameters') == asdict(p)
                and (output/(name+'_sources.npz')).exists()):
            print('Reusing completed refinement: '+name, flush=True)
            continue
        print('Smooth-wafer refinement: '+name, flush=True)
        points, weights, metrics = extract(build_mesh(p), p, True)
        np.savez_compressed(output/(name+'_sources.npz'), source_points=points,
                            source_polarization_weights=weights)
        report['cases'][name] = dict(parameters=asdict(p), solver=metrics)
        (output/'flat_refinement_native.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        del points, weights
        gc.collect()


def measure_refinements(output):
    from fullwave_prototype.far_field import radiation_amplitude
    from .channel_symmetry_audit import plot_pairs
    from .channel_symmetry_local import brightness
    native = json.loads((output/'flat_refinement_native.json').read_text(encoding='utf-8'))
    config = json.loads((Path(native['source'])/'hybrid_input.json').read_text(encoding='utf-8'))
    camera = DetectorParameters(**profile()['detector'])
    report, images = {}, []
    baseline = json.loads((output/'native_report.json').read_text(encoding='utf-8'))
    cases = dict(preview_order1=dict(parameters=baseline['parameters'], solver=baseline['flat']))
    cases.update(native['cases'])
    for name, case in cases.items():
        p = Parameters(**case['parameters'])
        settings = replace(optical_settings(config['imaging_settings'], p), padding=4)
        values, pictures = [], []
        filename = 'flat_sources.npz' if name == 'preview_order1' else name+'_sources.npz'
        with np.load(output/filename) as saved:
            for tilt in (45, -45):
                geometry = pupil_geometry(.5, tilt, settings.pupil_samples)
                amplitude = radiation_amplitude(saved['source_points'], saved['source_polarization_weights'],
                                                geometry[1][geometry[3]], p)
                metrics, image = channel_image(amplitude, geometry, p, settings, camera)
                values.append(metrics)
                pictures.append(image)
        errors = {key: abs(values[0][key]-values[1][key])/np.mean([v[key] for v in values])
                  for key in values[0]}
        errors['mirrored_image_relative_l2'] = relative_error(pictures[0]['flux'],
            reflect_centered_grid(pictures[1]['flux']))
        peaks = [brightness(image, camera)['central_peak_w_um2'] for image in pictures]
        errors['peak_symmetry_relative_difference'] = abs(peaks[0]-peaks[1])/np.mean(peaks)
        report[name] = dict(channels=values, peak_w_um2=peaks, errors=errors, solver=case['solver'])
        images.append(pictures)
    plot_pairs(output/'flat_refinement.png', images, list(report))
    (output/'flat_refinement_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('build/channel_symmetry_20261016'))
    parser.add_argument('--stage', choices=('solve', 'measure'), required=True)
    parser.add_argument('--fine-order2', action='store_true', help='Optional high-memory condensed refinement')
    args = parser.parse_args()
    if args.stage == 'solve':
        solve_refinements(args.output, args.fine_order2)
    else:
        measure_refinements(args.output)
