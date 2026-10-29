"""Assemble the periodic solver, corresponding source and original demo fields."""
import hashlib
import json
from pathlib import Path
import shutil


def sha256(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            value.update(block)
    return value.hexdigest()


def assemble_grating(app, root):
    from grating_scene.examples import EXAMPLES
    from grating_scene.camera import load_camera
    from grating_scene.pipeline import load_scene
    app, root = Path(app), Path(root)
    source = root/'build/s4_runtime'
    target = app/'s4'
    target.mkdir(parents=True, exist_ok=True)
    for name in ('S4.exe', 'openblas.dll', 'build_manifest.json'):
        shutil.copy2(source/'bin_stripes'/name, target/name)
    manifest = json.loads((target/'build_manifest.json').read_text(encoding='utf-8'))
    if sha256(target/'S4.exe') != manifest['executable_sha256']:
        raise ValueError('Bundled S4 executable does not match its build manifest')
    if sha256(target/'openblas.dll') != manifest['blas_sha256']:
        raise ValueError('Bundled OpenBLAS does not match its build manifest')
    shutil.copytree(source/'licenses', target/'licenses', dirs_exist_ok=True)
    corresponding = target/'source'
    corresponding.mkdir(exist_ok=True)
    shutil.copytree(source/'victorliu-S4-5b89848', corresponding/'victorliu-S4-5b89848',
        dirs_exist_ok=True, ignore=shutil.ignore_patterns('*.o', '*.a', '*.exe', '*.dll', '__pycache__'))
    for name in ('S4-5b89848cfd15.tar.gz', 'lua-5.2.4.tar.gz'):
        shutil.copy2(source/'downloads'/name, corresponding/name)
    for relative in ('tools/build_s4_runtime.py', 'tools/s4_native_extensions.py',
                     's4_prototype/patches/pre_api_gcc_compat.patch'):
        destination = corresponding/relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root/relative, destination)
    shutil.copy2(source/'stripe_extensions.patch', corresponding/'stripe_extensions.patch')
    shutil.copy2(root/'far_field_standalone/S4_SOURCE_README.md', corresponding/'README.md')
    examples = app/'examples/grating'
    examples.mkdir(parents=True, exist_ok=True)
    entries = []
    for row in EXAMPLES:
        original = root/row['source']
        (load_camera if row['kind'] == 'camera' else load_scene)(original)
        destination = examples/row['key']
        shutil.copytree(original, destination, dirs_exist_ok=True)
        (load_camera if row['kind'] == 'camera' else load_scene)(destination)
        entries.append(dict(key=row['key'], title=row['title'], kind=row['kind'],
            original_data_included=True, quantitative_hybrid_certified=False))
    (examples/'catalog.json').write_text(json.dumps(entries, indent=2), encoding='utf-8')
    for name in ('grating_contact_preview.json', 'grating_contact_camera_optics.json',
                 'grating_camera_gemi_area.json'):
        shutil.copy2(root/'examples'/name, app/'examples'/name)
    docs = app/'docs'
    docs.mkdir(exist_ok=True)
    for name in ('GRATING_SCENE_GUIDE_CN.md', 'GRATING_SCENE_DEVELOPMENT_CN.md',
                 'GRATING_SOLVER_ROADMAP_CN.md', 'S4_VALIDATION_CN.md'):
        shutil.copy2(root/'docs'/name, docs/name)
    development = docs/'GRATING_SCENE_DEVELOPMENT_CN.md'
    development.write_text(development.read_text(encoding='utf-8').replace(
        '../output/grating_scene_20261027/camera_gemi_area/camera.png',
        '../examples/grating/contact_camera_40x/camera.png'), encoding='utf-8')
    return dict(native_executable_sha256=sha256(target/'S4.exe'), examples=entries,
        quantitative_hybrid_certified=False, native_backend='S4, pinned stripes adapter',
        worker='fem/WaferScope_FEM.exe', corresponding_source='s4/source')


def copy_grating_verification(app, source):
    source, app = Path(source), Path(app)
    destination = app/'verification/grating_hybrid_validation'
    destination.mkdir(parents=True, exist_ok=True)
    for name in ('hybrid_summary.json', 'HYBRID_VALIDATION.md', 'tests.xml', 'ui_tests.xml',
                 'crosscheck_715.admission.json'):
        if (source/name).is_file():
            shutil.copy2(source/name, destination/name)
    for folder in source.glob('crosscheck_*'):
        if folder.is_dir() and (folder/'crosscheck.json').is_file():
            shutil.copytree(folder, destination/'raw'/folder.name, dirs_exist_ok=True)
            for suffix in ('.resources.json', '.admission.json'):
                path = folder.with_suffix(suffix)
                if path.is_file():
                    shutil.copy2(path, destination/'raw'/path.name)
    packaged = source/'packaged_r2'
    for name in ('grating_package_smoke.json', 'grating_near_ui.png', 'grating_camera_ui.png',
                 'grating_calculated_camera_ui.png'):
        if (packaged/name).is_file():
            target = destination/'packaged'/name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(packaged/name, target)
    for folder, name in (('scatmech_packaged', 'raw_export_smoke.json'),
                         ('relocated_packaged', 'grating_package_smoke.json')):
        path = source/folder/name
        if path.is_file():
            target = destination/folder/name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
