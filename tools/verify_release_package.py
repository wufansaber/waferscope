"""Check an extracted Windows release and collect current acceptance evidence."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET
import zipfile


ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def unpack(archive, destination):
    destination.mkdir(parents=True, exist_ok=False)
    with zipfile.ZipFile(archive) as saved:
        assert saved.testzip() is None, 'Corrupt release ZIP'
        for entry in saved.infolist():
            assert destination.resolve() in (destination/entry.filename).resolve().parents
        saved.extractall(destination)


def validate(app, built, evidence):
    for name in ('DFI_Simulation.exe', 'fem/DFI_FEM.exe'):
        assert digest(app/name) == digest(built/name), name
    for name in ('DFI_Simulation_中文操作手册.pdf', 'DFI_Simulation_功能与案例手册.pdf'):
        assert digest(app/'docs'/name) == digest(ROOT/'output/pdf'/name), name
    assert (app/'docs/far_field_standalone/manual_assets/fem_overview.png').is_file()
    assert (app/'docs/showcase_20261017/assets/case01_far_field_dense.png').is_file()
    tests = {}
    for name in ('nist2_tests.xml', 'native_tests.xml'):
        suites = ET.parse(evidence/name).getroot().iter('testsuite')
        counts = {key: 0 for key in ('tests', 'failures', 'errors', 'skipped')}
        for suite in suites:
            for key in counts:
                counts[key] += int(suite.get(key, 0))
        assert counts['tests'] > 0 and counts['failures'] == counts['errors'] == 0
        tests[name] = counts
    for name in ('packaged/package_smoke.json', 'packaged_camera/camera_ui_smoke.json'):
        assert json.loads((evidence/name).read_text(encoding='utf-8'))['status'] == 'passed'
    resources = json.loads((evidence/'cross_resources.json').read_text(encoding='utf-8'))
    assert resources['status'] == 'completed' and not resources['surviving_owned_pids']
    report = json.loads((evidence/'cross/imaging_report.json').read_text(encoding='utf-8'))
    assert report['scattering_parameters']['defect_shape'] == 'cross'
    assert all(ch['with']['image_power_w'] > 0 for ch in report['channels'])
    assert all(ch['detector']['gray_valid'] for ch in report['channels'])
    for name in ('near_field_xy.png', 'near_field_xz.png', 'near_field_yz.png',
                 'far_field.png', 'image_planes.png', 'detector.png', 'gray_acquisition.png'):
        assert (evidence/'cross'/name).is_file(), name
    assert len(list((evidence/'scatmech_dfi').glob('*.png'))) >= 2
    for name in ('manual_validation.json', 'showcase_validation.json'):
        assert json.loads((evidence/name).read_text(encoding='utf-8'))['status'] == 'passed'
    result = dict(status='passed', tests=tests, extracted_directory=str(app),
        gui_sha256=digest(app/'DFI_Simulation.exe'),
        fem_sha256=digest(app/'fem/DFI_FEM.exe'),
        cross_image_power_w=[ch['with']['image_power_w'] for ch in report['channels']],
        cross_gray_snr=[ch['detector']['gray_snr'] for ch in report['channels']],
        checks=['SCATMECH rendering', 'Frozen GUI to FEM calculation and cached re-imaging',
                'Frozen camera raw-DN SNR replay', 'Frozen finite cross near/far/image/camera',
                'Both PDF manuals match verified originals'],
        scope='Runtime and regression acceptance, not proof of FEM mesh/physical convergence.')
    (evidence/'package_validation.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    destination = built/'verification/release_20261019'
    for relative in ('package_validation.json', 'nist2_tests.xml', 'native_tests.xml',
                     'manual_validation.json', 'showcase_validation.json',
                     'packaged/package_smoke.json', 'packaged/packaged_gui.png',
                     'packaged_camera/camera_ui_smoke.json', 'packaged_camera/camera_metrics_1200.png',
                     'cross_resources.json', 'cross/imaging_report.json', 'cross/near_far_image.png'):
        target = destination/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(evidence/relative, target)
    print(json.dumps(result, indent=2))


def seal(archive, built):
    with zipfile.ZipFile(archive) as saved:
        assert saved.testzip() is None
        for entry in saved.infolist():
            if entry.is_dir():
                continue
            path = Path(entry.filename)
            assert path.parts[0] == 'DFI_Simulation'
            relative = Path(*path.parts[1:])
            assert '..' not in relative.parts
            assert hashlib.sha256(saved.read(entry)).hexdigest() == digest(built/relative), str(relative)
        assert 'DFI_Simulation/verification/release_20261019/package_validation.json' in saved.namelist()
    files = [archive, ROOT/'output/pdf/DFI_Simulation_中文操作手册.pdf',
             ROOT/'output/pdf/DFI_Simulation_功能与案例手册.pdf',
             ROOT/'output/pdf/DFI_Simulation_案例与宣传资料.zip',
             ROOT/'far_field_standalone/OPERATION_MANUAL_CN.md']
    checksums = archive.parent/'SHA256SUMS.txt'
    checksums.write_text(''.join(digest(path)+'  '+path.name+'\n' for path in files), encoding='utf-8')
    print(json.dumps(dict(status='passed', archive=str(archive), bytes=archive.stat().st_size,
                         sha256=digest(archive), checksums=str(checksums)), indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path)
    parser.add_argument('--extract-to', type=Path)
    parser.add_argument('--app', type=Path)
    parser.add_argument('--built', type=Path)
    parser.add_argument('--evidence', type=Path)
    parser.add_argument('--seal', type=Path)
    args = parser.parse_args()
    if args.seal:
        seal(args.seal, args.built)
    elif args.archive:
        unpack(args.archive, args.extract_to)
    else:
        validate(args.app, args.built, args.evidence)


if __name__ == '__main__':
    main()
