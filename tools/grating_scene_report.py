"""Rebuild the fixed contact-sphere development evidence from original records."""
import argparse
from pathlib import Path
import json
import xml.etree.ElementTree as ET

from grating_scene.camera import load_camera
from grating_scene.diagnostics import scene_diagnostics
from grating_scene.pipeline import digest, load_scene, write_json
from tools.grating_scene_validate import compare_scenes


def campaign(root):
    root = Path(root).resolve()
    base, _, _ = load_scene(root/'contact_preview')
    pairs = [('preview_x', 'contact_preview', 'native_controls/x_cutoff'),
        ('preview_y', 'contact_preview', 'native_controls/y_cutoff'),
        ('preview_slices', 'contact_preview', 'native_controls/slices'),
        ('preview_joint', 'contact_preview', 'native_controls/joint'),
        ('refined_x', 'native_controls/joint', 'native_controls_refined/x_cutoff'),
        ('refined_y', 'native_controls/joint', 'native_controls_refined/y_cutoff'),
        ('refined_slices', 'native_controls/joint', 'native_controls_refined/slices'),
        ('refined_joint', 'native_controls/joint', 'joint_profiled'),
        ('y32_to_y40', 'native_controls_refined/y_cutoff', 'y40_control')]
    comparisons = {name: compare_scenes(root/a, root/b) for name, a, b in pairs}
    scenes = {}
    for name in sorted({p for _, a, b in pairs for p in (a, b)}):
        config, report, arrays = load_scene(root/name)
        row = dict(report_sha256=digest(root/name/'scene_report.json'),
            dimensions=dict(mx=config.scene.harmonics_x, my=config.scene.harmonics_y, slices=config.scene.slices),
            collection=report['collection'], diagnostics=scene_diagnostics(arrays['far'], report['native'], config))
        resource = (root/name).with_suffix('.resources.json')
        row['resources'] = json.loads(resource.read_text(encoding='utf-8'))
        if row['resources']['status'] != 'completed' or row['resources']['surviving_owned_pids']:
            raise ValueError('Campaign contains an incomplete or unclean native task: '+name)
        scenes[name] = row
    camera, _ = load_camera(root/'camera_gemi_area')
    tests = {}
    for name in ('tests.xml', 'ui_tests.xml'):
        document = ET.parse(root/name).getroot()
        suites = [document] if document.tag == 'testsuite' else list(document.iter('testsuite'))
        totals = {k: sum(int(s.get(k, '0')) for s in suites) for k in ('tests', 'failures', 'errors', 'skipped')}
        if totals['failures'] or totals['errors']:
            raise ValueError('Regression is not passing: '+name)
        tests[name] = dict(totals, sha256=digest(root/name))
    smoke = json.loads((root/'qprocess_smoke_r3/qprocess_smoke.json').read_text(encoding='utf-8'))
    if smoke['status'] != 'passed':
        raise ValueError('Qt process smoke did not pass')
    return dict(status='periodic_workflow_implemented_not_isolated_certified',
        physical_scene=base.to_dict()['scene'], comparisons=comparisons, scenes=scenes,
        camera=dict(settings=camera['settings'], channels=camera['channels'],
            report_sha256=digest(root/'camera_gemi_area/camera_report.json'), optical_model_validated=False),
        regressions=tests, qprocess_smoke=smoke,
        limitations=['All new cases use gap=0; legacy gap=10 nm rounds are separate series.',
            'Local complex-field convergence and X mirror checks have not all passed.',
            'A periodic ideal-relay camera is implemented, not an isolated-defect real-lens certification.',
            'Plane wave, fixed uncalibrated material inputs, no haze/coatings in this scene.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Use a new summary filename')
    write_json(args.output, campaign(args.root))


if __name__ == '__main__':
    main()
