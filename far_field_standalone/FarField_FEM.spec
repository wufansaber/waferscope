# -*- mode: python ; coding: utf-8 -*-
from importlib.metadata import distribution
from pathlib import Path
import sys
import setuptools
from PyInstaller.utils.hooks import collect_all, collect_delvewheel_libs_directory

root = Path(SPECPATH).parent
datas, binaries, hiddenimports = [], [], ['backports.tarfile']
ffi = Path(sys.base_prefix)/'Library/bin/ffi.dll'
if ffi.is_file():
    binaries.append((str(ffi), '.'))
for package in ('ngsolve', 'netgen', 'pyngcore', 'ngsolve_openblas'):
    data, binary, hidden = collect_all(package, include_py_files=False)
    datas += [(src, dst) for src, dst in data if '/include/' not in src.replace('\\', '/')]
    binaries += binary
    hiddenimports += hidden
for package, library in (('ngsolve', 'ngsolve.libs'), ('netgen', 'netgen_mesher.libs')):
    datas, binaries = collect_delvewheel_libs_directory(package, library, datas=datas, binaries=binaries)
# OCCT's wheel RECORD uses ../../bin, not the frozen layout. Its metadata is
# intentionally omitted; native dependencies load from adjacent shared DLLs.
for file in distribution('netgen-occt').files:
    if str(file).lower().endswith('.dll'):
        binaries.append((str(file.locate()), '.'))
a = Analysis([str(Path(SPECPATH)/'fem_entry.py')],
    pathex=[str(root), str(Path(setuptools.__file__).parent/'_vendor')],
    binaries=binaries, datas=datas, hiddenimports=hiddenimports,
    excludes=['PyQt5', 'PySide6', 'tkinter', 'IPython', 'jupyter', 'numba', 'llvmlite', 'miepython'],
    noarchive=False)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name='DFI_FEM',
          console=True, debug=False, strip=False, upx=False)
coll = COLLECT(exe, a.binaries, a.datas, name='DFI_FEM', strip=False, upx=False)
