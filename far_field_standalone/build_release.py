"""Build and assemble the GUI and native solver as separate onedir executables."""
import argparse
from datetime import date
import hashlib
import importlib.metadata as metadata
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent


def export_licenses(destination):
    destination.mkdir(parents=True, exist_ok=True)
    names = ('numpy', 'scipy', 'matplotlib', 'pillow', 'pyparsing', 'contourpy',
        'cycler', 'kiwisolver', 'packaging', 'fonttools', 'python-dateutil', 'six',
        'setuptools', 'pyinstaller', 'PyQt5', 'PyQt5-Qt5', 'PyQt5-sip', 'pySCATMECH',
        'ngsolve', 'netgen-mesher', 'netgen-occt', 'ngsolve-openblas', 'psutil')
    entries = []
    for name in names:
        try:
            dist = metadata.distribution(name)
        except metadata.PackageNotFoundError:
            continue
        entries.append(dict(name=name, version=dist.version,
            homepage=dist.metadata.get('Home-page'),
            project_urls=dist.metadata.get_all('Project-URL', [])))
        for file in dist.files or []:
            if not any(word in file.name.lower() for word in ('license', 'copying', 'copyright', 'notice')):
                continue
            source = Path(file.locate())
            if source.is_file():
                relative = Path(*[part for part in file.parts if part not in ('.', '..')])
                target = destination/name/relative
                if destination.resolve() not in target.resolve().parents:
                    raise ValueError('License destination escaped output directory')
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
    (destination/'packages.json').write_text(json.dumps(entries, indent=2), encoding='utf-8')
    python_license = Path(sys.base_prefix)/'LICENSE.txt'
    if python_license.is_file():
        shutil.copy2(python_license, destination/'PYTHON_LICENSE.txt')


def copy_release_manuals(app, root=ROOT):
    """Keep Markdown image links and both reader-ready PDFs in the bundle."""
    for relative in ('far_field_standalone/manual_assets', 'docs/showcase_20261017'):
        source = root/relative
        if not source.is_dir():
            raise FileNotFoundError('Missing release documentation: '+str(source))
        destination = (app/'docs'/relative if relative.startswith('far_field_standalone/')
                       else app/relative)
        shutil.copytree(source, destination, dirs_exist_ok=True)
    for name in ('FarField_Simulation_中文操作手册.pdf', 'FarField_Simulation_功能与案例手册.pdf'):
        source = root/'output/pdf'/name
        if not source.is_file():
            raise FileNotFoundError('Missing release manual: '+str(source))
        destination = app/'docs'/name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


def build(args):
    output = HERE/'dist'/args.version
    app = output/'FarField_Simulation'
    engine = ROOT/'build/package_fem/dist/FarField_FEM'
    if not args.assemble_only:
        for python, spec, dist, work in (
            (args.gui_python, 'FarField_Simulation.spec', output, ROOT/'build/package_gui/work'),
            (args.fem_python, 'FarField_FEM.spec', engine.parent, ROOT/'build/package_fem/work')):
            subprocess.run([str(python), '-m', 'PyInstaller', str(HERE/spec), '--noconfirm',
                            '--distpath', str(dist), '--workpath', str(work)], cwd=str(HERE), check=True)
    if not (app/'FarField_Simulation.exe').is_file() or not (engine/'FarField_FEM.exe').is_file():
        raise FileNotFoundError('Build both executables before assembly')
    shutil.copytree(engine, app/'fem', dirs_exist_ok=True)
    for folder in ('far_field_standalone', 'fullwave_prototype'):
        target = app/'docs'/folder
        target.mkdir(parents=True, exist_ok=True)
        for source in (ROOT/folder).glob('*.md'):
            shutil.copy2(source, target/source.name)
    copy_release_manuals(app)
    shutil.copy2(HERE/'RELEASE_README_CN.md', app/'README_CN.md')
    shutil.copy2(HERE/'THIRD_PARTY_NOTICES.md', app/'THIRD_PARTY_NOTICES.md')
    for name in ('CAMERA_MODEL_EMVA1288_CN.md', 'DETECTION_SNR_UI_CN.md'):
        if (ROOT/name).is_file():
            shutil.copy2(ROOT/name, app/'docs'/name)
    for label, python in (('gui', args.gui_python), ('fem', args.fem_python)):
        subprocess.run([str(python), str(Path(__file__).resolve()), '--export-licenses',
                        str(app/'licenses'/label)], check=True)
    example = ROOT/'build/release_upgrade/reimage_final'
    if not (example/'imaging_report.json').is_file():
        example = ROOT/'build/near_pml_order/pipeline_outer'
    if (example/'imaging_report.json').is_file():
        target = app/'examples/validated_355nm_100nm'
        target.mkdir(parents=True, exist_ok=True)
        for name in ('imaging_report.json', 'near_field.png', 'far_field.png', 'image_planes.png',
                     'image_decomposition.png', 'near_far_image.png'):
            shutil.copy2(example/name, target/name)
    git = ['git', '-c', 'safe.directory='+ROOT.as_posix(), '-C', str(ROOT)]
    commit = subprocess.check_output(git+['rev-parse', 'HEAD'], text=True).strip()
    dirty = bool(subprocess.check_output(git+['status', '--porcelain'], text=True).strip())
    evidence = app/'verification'
    evidence.mkdir(exist_ok=True)
    audit_root = ROOT/'build/release_upgrade'
    for relative in ('sensitivity_summary.json', 'batch.json', 'beam/audit.json', 'beam_far/audit.json',
                     'reimage_final/audit.json', 'native_tests.json', 'packaged/package_smoke.json'):
        source = audit_root/relative
        if source.is_file():
            target = evidence/'historical_20260928'/relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
    if args.verification_root:
        for relative in ('campaign_summary.json', 'campaign_summary.png', 'pytest.txt',
                         'packaged/package_smoke.json', 'packaged/packaged_gui.png',
                         'packaged/whole_haze_smoke.json', 'packaged/haze_ensemble_smoke.json'):
            source = args.verification_root/relative
            if source.is_file():
                target = evidence/relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
        for name in ('summary_report.json', 'nist_haze_report.json', 'nist_haze_comparison.png',
                     'low_haze_comparison.json', 'low_haze_comparison.png', 'peak_decomposition.png',
                     'rough_mirror.png', 'ensemble_confidence.png', 'all_tests.xml',
                     'camera_ui_smoke.json', 'camera_metrics_1200.png', 'camera_metrics_640.png'):
            source = args.verification_root/name
            if source.is_file():
                target = evidence/'camera_20261016'/name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
        for name in ('camera_ui_smoke.json', 'camera_metrics_1200.png',
                     'camera_metrics_640.png', 'gray_acquisition_ui.png'):
            source = args.verification_root/'packaged_camera'/name
            if source.is_file():
                target = evidence/'camera_20261016'/'packaged'/name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
    (app/'build_manifest.json').write_text(json.dumps(dict(version=args.version,
        source_commit=commit, source_dirty=dirty,
        platform='Windows x64', format='onedir', fem_status='experimental',
        source='https://github.com/wufansaber/aptk',
        note='No local environments, credentials or computed field archives included.'), indent=2), encoding='utf-8')
    archive = Path(shutil.make_archive(str(output/('FarField_Simulation-'+args.version+'-win64')),
                                     'zip', root_dir=str(output), base_dir=app.name))
    digest = hashlib.sha256()
    with archive.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            digest.update(chunk)
    archive.with_suffix('.sha256').write_text(digest.hexdigest()+'  '+archive.name+'\n', encoding='ascii')
    print(json.dumps(dict(archive=str(archive), bytes=archive.stat().st_size,
                         sha256=digest.hexdigest()), indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gui-python', type=Path, default=Path('D:/anaconda/envs/NIST2/python.exe'))
    parser.add_argument('--fem-python', type=Path, default=ROOT/'build/ngsolve_env/Scripts/python.exe')
    parser.add_argument('--version', default=date.today().isoformat())
    parser.add_argument('--assemble-only', action='store_true')
    parser.add_argument('--verification-root', type=Path,
                        help='Current verification summaries to include, without field archives')
    parser.add_argument('--export-licenses', type=Path)
    args = parser.parse_args()
    if args.export_licenses:
        export_licenses(args.export_licenses)
    else:
        if not args.version or any(c not in '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ.-' for c in args.version) or args.version in ('.', '..'):
            parser.error('Use a simple version name without path separators')
        build(args)


if __name__ == '__main__':
    main()
