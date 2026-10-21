# -*- mode: python ; coding: utf-8 -*-
# WaferScope standalone GUI.

import os


project_root = os.path.abspath(os.path.join(SPECPATH, ".."))

a = Analysis(
    ["main.py"],
    pathex=[project_root],
    binaries=[],
    datas=[],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["pages", "main_window", "ngsolve", "netgen", "pyngcore", "ngsolve_openblas",
              "numba", "llvmlite", "miepython"],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="WaferScope",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="WaferScope",
)
