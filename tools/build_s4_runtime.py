"""Build a private Windows S4 Lua executable without touching Python envs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[1]
S4_COMMIT = '5b89848cfd15719ccaa04c853fe883ecf8d6ac15'
INPUTS = {
    'S4-5b89848cfd15.tar.gz': '136bd16e9ec0d16894e05a1ee15d1ee52fd1e507bfc84a0695da4b85c41ee2c3',
    'lua-5.2.4.tar.gz': 'b9e2e4aad6789b3b63a056d442f7b39f0ecfca3ae0f1fc0ae4e9614401b69f4b',
    'w64devkit-api.7z.exe': 'add1edf984ce4d473cfeb7c52e0ab38db73a855bba6c7f564d5e16648eb2896a',
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def extract_sources(archive, destination):
    destination = Path(destination).resolve()
    with tarfile.open(archive) as saved:
        members = saved.getmembers()
        for entry in members:
            (destination/entry.name).resolve().relative_to(destination)
            if entry.issym() or entry.islnk() or entry.isdev():
                raise ValueError('Source archive contains a link or special file')
        saved.extractall(destination, members=members, filter='data')


def build(root, blas, eigensolver='builtin', variant='baseline'):
    root, blas = Path(root).resolve(), Path(blas).resolve()
    for name, expected in INPUTS.items():
        if sha256(root/'downloads'/name) != expected:
            raise ValueError('Build input checksum mismatch: '+name)
    compiler = root/'w64devkit/bin'
    gcc, ar, make = [compiler/name for name in ('gcc.exe', 'ar.exe', 'make.exe')]
    if not all(path.is_file() for path in (gcc, ar, make, blas)):
        raise FileNotFoundError('Extract the checked toolchain and provide an OpenBLAS DLL')
    for name in ('S4-5b89848cfd15.tar.gz', 'lua-5.2.4.tar.gz'):
        extract_sources(root/'downloads'/name, root)
    source, lua = root/'victorliu-S4-5b89848', root/'lua-5.2.4'
    output = root/('bin_stripes' if variant == 'stripes' else 'bin')
    objects = root/'obj'
    output.mkdir(exist_ok=True)
    objects.mkdir(exist_ok=True)
    env = dict(os.environ, PATH=str(compiler)+os.pathsep+os.environ.get('PATH', ''),
               OPENBLAS_NUM_THREADS='2', OMP_NUM_THREADS='2')
    commands = []
    log = root/'build.log'
    log.write_text('', encoding='utf-8')

    def run(command, cwd=root):
        command = [str(value) for value in command]
        commands.append(command)
        with log.open('a', encoding='utf-8') as stream:
            stream.write('\n'+repr(command)+'\n')
            stream.flush()
            result = subprocess.run(command, cwd=cwd, env=env, stdout=stream,
                                    stderr=subprocess.STDOUT, timeout=240)
        if result.returncode:
            raise RuntimeError('Build failed; tail:\n'+'\n'.join(log.read_text(
                encoding='utf-8', errors='replace').splitlines()[-60:]))

    patch = ROOT/'s4_prototype/patches/pre_api_gcc_compat.patch'
    patch_directory = '--directory='+source.relative_to(ROOT).as_posix()
    run(['git', 'apply', '--unidiff-zero', patch_directory, '--check', patch], cwd=ROOT)
    run(['git', 'apply', '--unidiff-zero', patch_directory, patch], cwd=ROOT)
    extensions = None
    if variant == 'stripes':
        from tools.s4_native_extensions import apply_extensions
        generated_patch = root/'stripe_extensions.patch'
        header = apply_extensions(source, generated_patch)
        extensions = dict(header_sha256=sha256(header), patch_sha256=sha256(generated_patch),
            description='Project adapter: strict stripe Li material rule and explicit rectangular basis.')
    lua_objects = []
    for file in sorted((lua/'src').glob('*.c')):
        if file.name in ('lua.c', 'luac.c'):
            continue
        obj = objects/('lua_'+file.stem+'.o')
        run([gcc, '-c', '-O2', '-std=gnu99', '-DLUA_USE_WINDOWS', '-DLUA_COMPAT_ALL',
             '-I'+str(lua/'src'), file, '-o', obj])
        lua_objects.append(obj)
    lua_library = objects/'liblua52.a'
    run([ar, 'rcs', lua_library]+lua_objects)
    print('Lua static library built.', flush=True)
    shutil.copy2(blas, output/'openblas.dll')
    lapack = ' -DHAVE_LAPACK' if eigensolver == 'lapack' else ''
    cflags = '-DHAVE_BLAS'+lapack+' -O2 -std=gnu99 -fno-strict-aliasing'
    cxxflags = '-DHAVE_BLAS'+lapack+' -O2 -std=gnu++11 -fno-strict-aliasing'
    s4_objects = root/('objects_'+S4_COMMIT[:12]+'_'+eigensolver+'_'+variant)
    for folder in ('S4k', 'modules'):
        (s4_objects/folder).mkdir(parents=True, exist_ok=True)
    run([make, '-B', '-j2', 'OBJDIR='+s4_objects.as_posix(),
         'S4_LIBNAME='+(s4_objects/'libS4.a').as_posix(),
         'CC='+gcc.as_posix(), 'CXX='+(compiler/'g++.exe').as_posix(),
         'AR='+ar.as_posix(), 'CFLAGS='+cflags, 'CXXFLAGS='+cxxflags,
         'LUA_INC=-I'+(lua/'src').as_posix(),
         'LUA_LIB='+lua_library.as_posix()+' -static-libgcc -static-libstdc++',
         'LA_LIBS='+(output/'openblas.dll').as_posix(),
         (s4_objects/'S4').as_posix()], cwd=source)
    executable = s4_objects/'S4.exe'
    if not executable.is_file():
        executable = s4_objects/'S4'
    shutil.copy2(executable, output/'S4.exe')
    licenses = root/'licenses'
    licenses.mkdir(exist_ok=True)
    shutil.copy2(source/'COPYING', licenses/'S4_COPYING')
    shutil.copy2(source/'COPYRIGHT', licenses/'S4_COPYRIGHT')
    shutil.copy2(lua/'doc/readme.html', licenses/'Lua_readme.html')
    metadata = dict(backend='S4', source_commit=S4_COMMIT, frontend='Lua 5.2.4',
        variant=variant, native_extensions=extensions,
        eigensolver=eigensolver,
        source_url='https://github.com/victorliu/S4', inputs=INPUTS,
        compatibility_patch_sha256=sha256(patch),
        source_selection='Official pre-C-API-rewrite commit; thickness-pointer validation fix only.',
        compiler=subprocess.check_output([str(gcc), '--version'], env=env, text=True).splitlines()[0],
        executable_sha256=sha256(output/'S4.exe'), blas_sha256=sha256(output/'openblas.dll'),
        commands=commands, status='built_not_yet_physically_validated')
    if extensions is not None:
        metadata['source_selection'] = 'Pinned official S4 plus documented project numerical adapter.'
    encoded = json.dumps(metadata, indent=2)
    (output/'build_manifest.json').write_text(encoded, encoding='utf-8')
    if variant == 'baseline':
        (root/'build_manifest.json').write_text(encoded, encoding='utf-8')
    print(json.dumps(dict(executable=str(output/'S4.exe'), status=metadata['status']), indent=2))
    return metadata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT/'build/s4_runtime')
    parser.add_argument('--blas', type=Path,
        default=ROOT/'build/ngsolve_env/Lib/site-packages/ngsolve_openblas/openblas.dll')
    parser.add_argument('--eigensolver', choices=('builtin', 'lapack'), default='builtin')
    parser.add_argument('--variant', choices=('baseline', 'stripes'), default='baseline')
    args = parser.parse_args()
    build(args.root, args.blas, args.eigensolver, args.variant)
