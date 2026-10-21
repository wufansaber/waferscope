"""Create a source migration ZIP from a committed Git tree, not a dirty checkout."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import subprocess
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'DFI_Simulation_Migration_20261021'
EXCLUDED = {'.git', '.codex', '.agents', '.idea', '.vscode', '__pycache__', 'build',
            'dist', '.pytest_cache', '.mplconfig', '.venv', 'venv', 'env', 'offline_wheels'}


def allowed(name):
    path = PurePosixPath(name)
    return (not path.is_absolute() and '..' not in path.parts and '\\' not in name
            and not EXCLUDED.intersection(path.parts) and path.name != '.env'
            and not path.name.startswith('.env.')
            and path.suffix.lower() not in ('.pem', '.pfx', '.key', '.pyc'))


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def git(*arguments):
    return subprocess.check_output(['git', '-c', 'safe.directory='+ROOT.as_posix(),
        '-C', str(ROOT), *arguments])


def environment_snapshot(python):
    names = ['numpy', 'scipy', 'matplotlib', 'PyQt5', 'PyQt5-Qt5', 'PyQt5-sip',
        'pySCATMECH', 'ngsolve', 'netgen-mesher', 'netgen-occt', 'ngsolve-openblas',
        'psutil', 'pytest', 'miepython', 'grcwa', 'autograd', 'pyinstaller']
    probe = ("import importlib.metadata as m,json,platform,sys; "
        "names="+repr(names)+"; installed={d.metadata['Name'].lower():d.version for d in m.distributions()}; "
        "print(json.dumps(dict(python=sys.version,platform=platform.platform(),"
        "packages={n:installed.get(n.lower()) for n in names})))")
    return json.loads(subprocess.check_output([str(python), '-c', probe], text=True, timeout=60))


def write_bundle(source_zip, target, metadata, environments, omissions):
    records, excluded = [], []
    with zipfile.ZipFile(source_zip) as source, zipfile.ZipFile(target, 'x', zipfile.ZIP_DEFLATED) as out:
        for entry in source.infolist():
            if entry.is_dir():
                continue
            if not allowed(entry.filename):
                excluded.append(entry.filename)
                continue
            data = source.read(entry)
            name = PREFIX+'/'+entry.filename
            out.writestr(name, data)
            records.append(dict(path=entry.filename, bytes=len(data), sha256=sha256(data)))
        data = json.dumps(environments, indent=2).encode('utf-8')
        out.writestr(PREFIX+'/ENVIRONMENT_SNAPSHOTS.json', data)
        records.append(dict(path='ENVIRONMENT_SNAPSHOTS.json', bytes=len(data), sha256=sha256(data)))
        manifest = dict(**metadata, files=records, excluded_tracked_paths=excluded,
            omitted_local_outputs=omissions,
            note='Source-only migration; no EXE, portable venv, credentials or Git history. '
                 'The manifest does not hash itself. Omitted outputs remain on the original PC.')
        out.writestr(PREFIX+'/MIGRATION_MANIFEST.json', json.dumps(manifest, indent=2).encode('utf-8'))
    verify_bundle(target)
    return manifest


def verify_bundle(path):
    with zipfile.ZipFile(path) as archive:
        if archive.testzip() is not None:
            raise ValueError('ZIP CRC check failed')
        manifest = json.loads(archive.read(PREFIX+'/MIGRATION_MANIFEST.json'))
        expected = {PREFIX+'/MIGRATION_MANIFEST.json'}
        for record in manifest['files']:
            if not allowed(record['path']):
                raise ValueError('Excluded path in migration package')
            name = PREFIX+'/'+record['path']
            data = archive.read(name)
            if len(data) != record['bytes'] or sha256(data) != record['sha256']:
                raise ValueError('Migration content hash mismatch: '+name)
            expected.add(name)
        if set(archive.namelist()) != expected or len(archive.namelist()) != len(expected):
            raise ValueError('Unlisted or duplicated ZIP entries')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--ref', default='HEAD')
    parser.add_argument('--gui-python', type=Path)
    parser.add_argument('--fem-python', type=Path)
    parser.add_argument('--rcwa-python', type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    if ROOT not in output.parents:
        parser.error('Output must be inside this workspace')
    output.mkdir(parents=True, exist_ok=True)
    commit = git('rev-parse', '--verify', args.ref+'^{commit}').decode().strip()
    tree = git('rev-parse', commit+'^{tree}').decode().strip()
    environments = {label: environment_snapshot(python) for label, python in (
        ('gui_reference', args.gui_python), ('fem', args.fem_python), ('rcwa', args.rcwa_python)) if python}
    target = output/(PREFIX+'.zip')
    with tempfile.TemporaryDirectory(prefix='migration_', dir=output) as temporary:
        source = Path(temporary)/'source.zip'
        git('archive', '--format=zip', '--output='+str(source), commit)
        with zipfile.ZipFile(source) as archive:
            included = set(archive.namelist())
        omissions = [dict(path=file.relative_to(ROOT).as_posix(), bytes=file.stat().st_size,
            sha256=sha256(file.read_bytes())) for file in sorted((ROOT/'output').rglob('*'))
            if file.is_file() and file.relative_to(ROOT).as_posix() not in included]
        manifest = write_bundle(source, target, dict(source_commit=commit, source_tree=tree,
            repository='https://github.com/wufansaber/aptk', date='2026-10-21'), environments, omissions)
    digest = sha256(target.read_bytes())
    target.with_suffix('.sha256').write_text(digest+'  '+target.name+'\n', encoding='ascii')
    (output/'MIGRATION_MANIFEST.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    (output/'ENVIRONMENT_SNAPSHOTS.json').write_text(json.dumps(environments, indent=2), encoding='utf-8')
    print(json.dumps(dict(archive=str(target), bytes=target.stat().st_size, sha256=digest,
        source_commit=commit, files=len(manifest['files']), omitted_outputs=len(omissions)), indent=2))


if __name__ == '__main__':
    main()
