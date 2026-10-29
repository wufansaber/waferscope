"""Archive the existing Windows S4 runtime and corresponding source; no solves."""
import argparse
import json
from pathlib import Path
import zipfile

from tools.build_migration_bundle import ROOT, sha256_file

NAME = 'WaferScope_S4_Modified_20261029'


def build(output):
    runtime = ROOT/'build/s4_runtime'
    native = runtime/'bin_stripes'
    manifest = json.loads((native/'build_manifest.json').read_text(encoding='utf-8'))
    if manifest.get('variant') != 'stripes':
        raise ValueError('The explicit rectangular basis runtime is required')
    for name, key in [('S4.exe', 'executable_sha256'), ('openblas.dll', 'blas_sha256')]:
        if sha256_file(native/name) != manifest[key]:
            raise ValueError('Native runtime checksum mismatch: '+name)
    for name, expected in manifest['inputs'].items():
        if name == 'w64devkit-api.7z.exe':
            continue  # The compiler is downloaded separately and is not repackaged.
        if sha256_file(runtime/'downloads'/name) != expected:
            raise ValueError('Build input checksum mismatch: '+name)
    patch = ROOT/'s4_prototype/patches/pre_api_gcc_compat.patch'
    if sha256_file(patch) != manifest['compatibility_patch_sha256']:
        raise ValueError('Compatibility patch does not match the native build')
    if sha256_file(runtime/'stripe_extensions.patch') != manifest['native_extensions']['patch_sha256']:
        raise ValueError('Native adapter patch does not match the native build')
    source = runtime/'victorliu-S4-5b89848'
    if sha256_file(source/'S4/waferscope_fmm.h') != manifest['native_extensions']['header_sha256']:
        raise ValueError('Native adapter source does not match the native build')
    files = [native/name for name in ('S4.exe', 'build_manifest.json')]
    files += sorted(p for p in (runtime/'licenses').rglob('*') if p.is_file())
    files += [runtime/'downloads'/name for name in manifest['inputs'] if name != 'w64devkit-api.7z.exe']
    files += [runtime/'stripe_extensions.patch']
    files += sorted(p for p in source.rglob('*') if p.is_file()
                    and p.suffix.lower() not in ('.o', '.a', '.exe', '.dll', '.pyc')
                    and '__pycache__' not in p.parts)
    files += [ROOT/name for name in ('tools/build_s4_runtime.py', 'tools/s4_native_extensions.py',
        's4_prototype/native/waferscope_fmm.h', 's4_prototype/patches/pre_api_gcc_compat.patch',
        'far_field_standalone/S4_SOURCE_README.md', 'docs/PROJECT_MIGRATION_CN.md')]
    for license_name in ('OpenBLAS_LICENSE.txt', 'LAPACK_LICENSE.txt', 'GCC_RUNTIME_EXCEPTION.txt'):
        if not (runtime/'licenses'/license_name).is_file():
            raise FileNotFoundError('Add upstream license before distribution: '+license_name)
    output = Path(output).resolve()
    output.relative_to(ROOT)
    output.mkdir(parents=True, exist_ok=True)
    target = output/(NAME+'.zip')
    records = []
    with zipfile.ZipFile(target, 'x', zipfile.ZIP_DEFLATED) as archive:
        for path in files:
            relative = path.relative_to(ROOT).as_posix()
            archive.write(path, NAME+'/'+relative)
            records.append(dict(path=relative, bytes=path.stat().st_size, sha256=sha256_file(path)))
        inventory = dict(format='waferscope-s4-runtime-migration', date='2026-10-29',
            platform='Windows x64', source_commit=manifest['source_commit'], files=records,
            note='Modified S4 binary and matching corresponding source. OpenBLAS and compiler are '
                 'downloaded separately from official sources; see PROJECT_MIGRATION_CN.md. '
                 'No Python environment or simulation data. Extract into checkout root.')
        archive.writestr(NAME+'/S4_MIGRATION_MANIFEST.json', json.dumps(inventory, indent=2))
    # Stream ZIP CRC verification and compare the archived bytes with every listed hash.
    verify(target)
    digest = sha256_file(target)
    target.with_suffix('.sha256').write_text(digest+'  '+target.name+'\n', encoding='ascii')
    (output/'S4_MIGRATION_MANIFEST.json').write_text(json.dumps(inventory, indent=2), encoding='utf-8')
    print(json.dumps(dict(archive=str(target), bytes=target.stat().st_size,
        sha256=digest, files=len(records), native_executable_sha256=manifest['executable_sha256']), indent=2))


def verify(path):
    import hashlib
    with zipfile.ZipFile(path) as archive:
        if archive.testzip() is not None:
            raise ValueError('Runtime ZIP CRC check failed')
        inventory = json.loads(archive.read(NAME+'/S4_MIGRATION_MANIFEST.json'))
        expected = {NAME+'/S4_MIGRATION_MANIFEST.json'}
        for row in inventory['files']:
            relative = Path(row['path'])
            if relative.is_absolute() or '..' in relative.parts:
                raise ValueError('Unsafe runtime archive path')
            name = NAME+'/'+row['path']
            expected.add(name)
            value, size = hashlib.sha256(), 0
            with archive.open(name) as stream:
                for block in iter(lambda: stream.read(1024*1024), b''):
                    value.update(block)
                    size += len(block)
            if size != row['bytes'] or value.hexdigest() != row['sha256']:
                raise ValueError('Runtime archive hash mismatch: '+name)
        if len(archive.namelist()) != len(expected) or set(archive.namelist()) != expected:
            raise ValueError('Unlisted or duplicate runtime archive entries')
    return inventory


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    build(parser.parse_args().output)
