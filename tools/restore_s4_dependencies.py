"""Download the pinned OpenBLAS wheel and restore only the DLL required by S4."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile

from grating_scene.runtime import default_s4


def restore(wheel, target):
    target = Path(target)
    manifest = json.loads((target/'build_manifest.json').read_text(encoding='utf-8'))
    expected = manifest['blas_sha256']
    destination = target/'openblas.dll'
    if destination.exists():
        if hashlib.sha256(destination.read_bytes()).hexdigest() == expected:
            print('Matching OpenBLAS DLL is already present.')
            return
        raise FileExistsError('Existing OpenBLAS has a different hash; preserve it before retrying')
    with zipfile.ZipFile(wheel) as archive:
        members = [n for n in archive.namelist() if n == 'ngsolve_openblas/openblas.dll']
        if len(members) != 1:
            raise ValueError('A Windows x64 ngsolve-openblas 0.3.33 wheel is required')
        with tempfile.TemporaryDirectory(prefix='blas_', dir=target) as temporary:
            candidate = Path(temporary)/'openblas.dll'
            with archive.open(members[0]) as source, candidate.open('wb') as output:
                shutil.copyfileobj(source, output, length=1024*1024)
            if hashlib.sha256(candidate.read_bytes()).hexdigest() != expected:
                raise ValueError('Downloaded DLL does not match the pinned S4 build')
            candidate.rename(destination)
    print('OpenBLAS restored; SHA-256: '+expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wheel', type=Path, help='Previously downloaded Windows x64 wheel')
    parser.add_argument('--download', action='store_true', help='Download the pinned wheel from PyPI')
    parser.add_argument('--target', type=Path, default=default_s4().parent)
    args = parser.parse_args()
    if bool(args.wheel) == bool(args.download):
        parser.error('Choose --download or --wheel')
    if not (args.target/'build_manifest.json').is_file():
        parser.error('Extract the modified S4 package into this checkout first')
    if args.download:
        with tempfile.TemporaryDirectory(prefix='s4_wheel_') as temporary:
            subprocess.run([sys.executable, '-m', 'pip', 'download', '--index-url',
                'https://pypi.org/simple', '--only-binary=:all:', '--no-deps',
                '--platform', 'win_amd64', '--python-version', '311',
                '--implementation', 'cp', '--abi', 'cp311', '--dest', temporary,
                'ngsolve-openblas==0.3.33'], check=True)
            wheels = list(Path(temporary).glob('*.whl'))
            if len(wheels) != 1:
                raise ValueError('Expected exactly one dependency wheel')
            restore(wheels[0], args.target)
    else:
        restore(args.wheel, args.target)


if __name__ == '__main__':
    main()
