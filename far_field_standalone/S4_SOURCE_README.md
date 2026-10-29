# Corresponding S4 Source

The adjacent `S4.exe` is a project-specific build of Stanford S4 commit
`5b89848cfd15719ccaa04c853fe883ecf8d6ac15`, with a GCC compatibility fix and
the explicit rectangular-basis / stripe Li-rule adapter. It uses the built-in
RNP eigensolver and a separately distributed OpenBLAS DLL. It is not an
unmodified upstream binary. The GPL v2-or-later notice and full license are
in `../licenses/` and the source tree.

This directory includes the modified S4 source tree, pristine upstream S4 and
Lua 5.2.4 archives, both patches, the generator, and the Python build recipe.
`../build_manifest.json` records input hashes, compiler version, all compilation
commands, and output hashes. Lua's copyright/license is in its archive and
`../licenses/Lua_readme.html`. OpenBLAS notices are included in the distribution's
`licenses/fem` inventory from the supplying `ngsolve-openblas` wheel.

To reproduce the project's Windows build, use the recipe in
`tools/build_s4_runtime.py` from the WaferScope repository. Its `--root` directory
must contain `downloads/` with the checked S4/Lua/toolchain archives and an
extracted `w64devkit/bin` toolchain; the toolchain archive is not part of this
end-user bundle. Obtain the checked compiler distribution using the recipe's
recorded input hash. Use `--variant stripes --eigensolver builtin --blas PATH`
with the matching OpenBLAS DLL. The recipe extracts pristine sources, applies
the included compatibility patch, generates the stripe adapter, builds Lua and
S4, and writes the runtime manifest. The exact compiler commands are retained
even when checkout paths differ. Compiler/toolchain dependencies are build-time
tools, not required to run the bundled executable.

Project source and recipes: https://github.com/wufansaber/waferscope

Upstream source: https://github.com/victorliu/S4

S4 is distributed without warranty. Supplying source and license information
does not imply the numerical model is converged, nor that this modified build
has upstream certification.
