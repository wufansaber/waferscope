"""Resource-bounded 2 um cross FEM comparison with a recorded saved baseline."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from fullwave_prototype.imaging_options import DETECTOR_OPTIONS


def command(source, output, scale, solver, particle_mesh=.03):
    r = json.loads((source/'imaging_report.json').read_text(encoding='utf-8'))
    p, s = r['scattering_parameters'], r['settings']
    values = dict(output=str(output), **{'defect-shape': 'cross', 'defect-length-nm': 2000,
        'defect-width-nm': 500, 'defect-height-nm': p['defect_height_nm'],
        'defect-rotation-deg': p['defect_rotation_deg'], 'gap-nm': p['gap_nm'],
        'domain-radius-um': 1.25, 'mesh-scale': scale, 'particle-mesh-um': particle_mesh,
        'accuracy': 'preview', 'order': 1, 'threads': 2, 'direct-solver': solver,
        'rms-nm': 0, 'roughness-model': 'gaussian', 'rough-grid-points': 9,
        'beam-profile': 'plane', 'illumination-fit': 'manual', 'pupil-samples': 64,
        'padding': 8, 'archive-mode': 'compact', 'pml-strength': p['pml_strength'],
        'pml-power': p['pml_power'], 'pml-mesh-scale': p['pml_mesh_scale'],
        'magnification': s['magnification'], 'objective-focal-mm': s['focal_length_mm'],
        'transmission': s['transmission'], 'defocus-um': s['image_defocus_um'],
        'analyzer': s['analyzer'], 'camera-profile': s.get('camera_profile', 'custom'),
        'gray-snr-pairs': 100, 'gray-snr-seed': s.get('gray_snr_seed', 20261016)})
    for key in ('wavelength_nm', 'particle_n', 'particle_k', 'substrate_n', 'substrate_k',
                'incidence_deg', 'polarization'):
        values[key.replace('_', '-')] = p[key]
    for key in ('power_w', 'length_um', 'width_um'):
        values['illumination-'+key.replace('_', '-')] = s['illumination_'+key]
    for key in DETECTOR_OPTIONS:
        values['detector-'+key] = s['detector'][key.replace('-', '_')]
    result = ['-m', 'fullwave_prototype.imaging_run']
    for key, value in values.items():
        result.extend(('--'+key, str(value)))
    result += ['--na']+list(map(str, s['numerical_apertures']))
    result += ['--tilts']+list(map(str, s['tilt_degrees']))
    if s['detector']['enabled']:
        result.append('--detector')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, default=Path('build/fullwave_ui/20261017_164144_979822'))
    parser.add_argument('--output', type=Path, default=Path('build/cross_2um_20261017'))
    parser.add_argument('--mesh-scale', type=float, default=3.)
    parser.add_argument('--particle-mesh-um', type=float, default=.03)
    parser.add_argument('--solver', default='sparsecholesky', choices=('umfpack', 'sparsecholesky'))
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable]+command(args.baseline, output, args.mesh_scale, args.solver, args.particle_mesh_um)
    (output/'comparison_request.json').write_text(json.dumps(dict(baseline=str(args.baseline.resolve()),
        command=cmd, note='Larger physical defect; coarser preview FEM, not a mesh-converged accuracy claim'),
        indent=2), encoding='utf-8')
    # Capture the worker's native messages while retaining the normal process-tree guard.
    worker = [sys.executable, '-m', 'fullwave_prototype.guarded_run', '--report', str(output/'resources.json'),
        '--max-memory-gib', '10', '--timeout-seconds', '900', 'imaging_run']+cmd[3:]
    with (output/'solve.log').open('w', encoding='utf-8') as stream:
        completed = subprocess.run(worker, stdout=stream, stderr=subprocess.STDOUT, timeout=960)
    if completed.returncode:
        raise RuntimeError('FEM comparison failed; inspect '+str(output/'resources.json'))
    report = json.loads((output/'imaging_report.json').read_text(encoding='utf-8'))
    print(json.dumps(dict(status='completed', output=str(output),
        power_w=[ch['with']['image_power_w'] for ch in report['channels']],
        gray_snr=[ch['detector']['gray_snr'] for ch in report['channels']],
        gray_valid=[ch['detector']['gray_valid'] for ch in report['channels']]), indent=2))


if __name__ == '__main__':
    main()
