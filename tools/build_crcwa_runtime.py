"""Private Windows candidate build; no solver acceptance or UI replacement implied."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess


COMMIT = '446a3a6ea5697db0bbb94af06d5616144b508c5a'
INPUTS = {'crcwa-source.tar.gz': '4adc66986479e71b5cd5067e19eb50cfe89016e8855c40434bab61a9cd234239',
          'fftw-3.3.5-dll64.zip': 'cfd88dc0e8d7001115ea79e069a2c695d52c8947f5b4f3b7ac54a192756f439f'}
SOURCES = ('block_matrix', 'main', 'matrixDev', 'parsefile', 'structure', 'section',
    'commands_structure', 'commands_section', 'commands_propagation', 'commands_other',
    'finterface', 'parsercl', 'matrixdevsym', 'parallelize_POSIX', 'toeplitz')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('build/crcwa_runtime_20261027'))
    args = parser.parse_args()
    root = args.root.resolve()
    workspace = Path(__file__).resolve().parents[1]
    compiler = workspace/'build/s4_runtime/w64devkit/bin/g++.exe'
    blas = workspace/'build/s4_runtime/bin_stripes/openblas.dll'
    source = root/('cRCWA-'+COMMIT)/'src'
    fftw = root/'fftw'
    for name, expected in INPUTS.items():
        if digest(root/name) != expected:
            raise ValueError('Archive checksum mismatch: '+name)
    directory = '--directory='+source.parent.relative_to(workspace).as_posix()
    patches = [workspace/'rcwa_prototype/patches'/name for name in
        ('crcwa_mingw_compat.patch', 'crcwa_binary_script_offsets.patch')]
    for patch in patches:
        check = subprocess.run(['git', 'apply', '--unidiff-zero', directory, '--check', str(patch)],
                                cwd=workspace, capture_output=True)
        if check.returncode == 0:
            subprocess.run(['git', 'apply', '--unidiff-zero', directory, str(patch)], cwd=workspace, check=True)
        else:
            subprocess.run(['git', 'apply', '--unidiff-zero', directory, '--reverse', '--check', str(patch)],
                           cwd=workspace, check=True)
    output, objects = root/'bin', root/'objects'
    output.mkdir(exist_ok=True)
    objects.mkdir(exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    log = root/'build.log'
    commands = []
    manifest = dict(source_url='https://github.com/cRCWA/cRCWA', source_commit=COMMIT,
        input_sha256=INPUTS, checksum_scope='Retrieved from official HTTPS endpoints and frozen; not a signature audit.',
        source_modifications=[dict(patch=str(patch), sha256=digest(patch)) for patch in patches],
        commands=commands, status='building',
        license='GPL-3.0-or-later: candidate kept separate from the application package.',
        acceptance='Build and small benchmarks only; no AFMM isolated-particle acceptance.')
    try:
        with log.open('w', encoding='utf-8') as stream:
            def run(command):
                command = [str(a) for a in command]
                commands.append(command)
                stream.write(repr(command)+'\n')
                stream.flush()
                subprocess.run(command, cwd=source, env=env, stdout=stream,
                               stderr=subprocess.STDOUT, timeout=180, check=True)
            paths = []
            for name in SOURCES:
                obj = objects/(name+'.o')
                run([compiler, '-std=gnu++11', '-O2', '-g0', '-D', 'LAPACK_ADD_FINAL_UNDERSCORE',
                    '-D', 'BLAS_ADD_FINAL_UNDERSCORE', '-D', 'USE_UNNAMED',
                    '-I'+str(fftw), '-c', source/(name+'.cpp'), '-o', obj])
                paths.append(obj)
                print('Compiled '+name, flush=True)
            executable = output/'crcwa.exe'
            run([compiler]+paths+[blas, fftw/'libfftw3-3.dll', '-pthread',
                '-static-libgcc', '-static-libstdc++', '-o', executable])
        shutil.copy2(blas, output/'openblas.dll')
        shutil.copy2(fftw/'libfftw3-3.dll', output/'libfftw3-3.dll')
        manifest.update(status='built_not_physically_validated', executable=str(executable),
            executable_sha256=digest(executable), blas_sha256=digest(blas),
            fftw_dll_sha256=digest(fftw/'libfftw3-3.dll'),
            source_files_sha256={name: digest(source/name) for name in
                                 [n+'.cpp' for n in SOURCES]+['compat.h']})
    except Exception as error:
        manifest.update(status='build_failed', error=str(error))
        raise
    finally:
        (root/'build_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == '__main__':
    main()
