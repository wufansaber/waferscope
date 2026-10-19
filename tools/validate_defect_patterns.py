"""Native CAD checks and bounded near/far/image/camera solves for finite patterns."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from fullwave_prototype.defects import cad_solid, metadata
from fullwave_prototype.parameters import Parameters


def cad_checks():
    import ngsolve as ng
    from fullwave_prototype.solver import build_mesh
    result = {}
    for shape in ('line_x', 'line_y', 'cross'):
        p = Parameters(defect_shape=shape, defect_rotation_deg=17, gap_nm=10)
        expected = metadata(p)['expected_volume_um3']
        actual = cad_solid(p, 1000).mass/1e9
        if not np.isclose(actual, expected, rtol=1e-12, atol=0):
            raise RuntimeError('Incorrect CAD union volume: '+shape)
        result[shape] = dict(expected_volume_um3=expected, cad_volume_um3=actual)
    p = replace(p, defect_rotation_deg=0, order=1, particle_mesh_um=.02,
                air_mesh_um=.12, substrate_mesh_um=.06)
    mesh = build_mesh(p)
    volume = float(ng.Integrate(1, mesh, definedon=mesh.Materials('particle')))
    if not np.isclose(volume, metadata(p)['expected_volume_um3'], rtol=1e-9, atol=0):
        raise RuntimeError('Incorrect cross mesh material volume')
    result['cross_mesh'] = dict(material_volume_um3=volume, elements=mesh.ne)
    return result


def run_case(root, shape, rough=False):
    name = shape+('_rough' if rough else '')
    folder = root/name
    command = [sys.executable, '-m', 'fullwave_prototype.guarded_run',
        '--report', str(root/(name+'_resources.json')), '--max-memory-gib', '8',
        '--timeout-seconds', '300', 'imaging_run', '--output', str(folder),
        '--defect-shape', shape, '--defect-length-nm', '200', '--defect-width-nm', '50',
        '--defect-height-nm', '50', '--wavelength-nm', '266', '--gap-nm', '10',
        '--rms-nm', '.2' if rough else '0', '--roughness-model', 'gaussian',
        '--rough-grid-points', '5' if rough else '9', '--correlation-length-um', '.1' if rough else '.05',
        '--accuracy', 'preview', '--order', '1',
        '--mesh-scale', '1.5', '--particle-mesh-um', '.015', '--threads', '2',
        '--direct-solver', 'umfpack', '--beam-profile', 'plane', '--illumination-fit', 'manual',
        '--illumination-power-w', '1', '--illumination-length-um', '2', '--illumination-width-um', '2',
        '--pupil-samples', '32', '--padding', '4', '--archive-mode', 'compact',
        '--transmission', '.000001', '--detector', '--detector-pixel-pitch-um', '5',
        '--detector-pixels-x', '20', '--detector-pixels-y', '20', '--detector-exposure-ms', '.01',
        '--detector-gain-e-dn', '16', '--detector-offset-dn', '64', '--detector-bit-depth', '10',
        '--detector-full-well-e', '16000', '--detector-read-noise-e', '12.1',
        '--detector-quantum-efficiency', '.639', '--gray-snr-pairs', '100']
    with (root/(name+'.log')).open('w', encoding='utf-8') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=340)
    scattering = json.loads((folder/'scattering/report.json').read_text(encoding='utf-8'))
    imaging = json.loads((folder/'imaging_report.json').read_text(encoding='utf-8'))
    expected = metadata(scattering['parameters'])['expected_volume_um3']
    actual = scattering['metrics']['with']['particle_volume_um3']
    if not np.isclose(actual, expected, rtol=1e-8, atol=0):
        raise RuntimeError('FEM material volume mismatch: '+name)
    for filename in ('near_field_xz.png', 'near_field_xy.png', 'near_field_yz.png',
                     'far_field.png', 'image_planes.png', 'detector.png', 'gray_acquisition.png'):
        if not (folder/filename).is_file():
            raise RuntimeError('Missing result: '+filename)
    channels = [ch['with']['image_power_w'] for ch in imaging['channels']]
    if not np.isfinite(channels).all() or min(channels) <= 0:
        raise RuntimeError('Invalid image powers')
    if not all(ch['detector']['gray_valid'] for ch in imaging['channels']):
        raise RuntimeError('Gray acquisition is clipped or invalid')
    return dict(status='passed', shape=shape, roughness_rms_nm=.2 if rough else 0,
        correlation_length_um=scattering['roughness']['correlation_length_um'],
        rough_grid_points=scattering['roughness']['grid_points'],
        expected_volume_um3=expected, fem_volume_um3=actual, image_power_w=channels,
        gray_snr=[ch['detector']['gray_snr'] for ch in imaging['channels']],
        solver_residual=scattering['metrics']['with']['algebraic_residual'],
        dofs=scattering['metrics']['with']['dofs'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('build/defect_patterns_20261017'))
    parser.add_argument('--include-rough', action='store_true')
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    report = dict(cad=cad_checks(), cases=[],
        scope='Preview solver / geometry integration checks; not a grid or sharp-corner convergence proof')
    for shape, rough in [('line_x', False), ('line_y', False), ('cross', False)]+(
            [('cross', True)] if args.include_rough else []):
        report['cases'].append(run_case(root, shape, rough))
        (root/'pattern_validation.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        print(json.dumps(report['cases'][-1]), flush=True)


if __name__ == '__main__':
    main()
