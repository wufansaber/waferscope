"""Call the private S4 Lua executable and retain all raw calculation records."""
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np

from .lua_scene import script
from .output import parse, reconstruct

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXE = ROOT/'build/s4_runtime/bin/S4.exe'


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def solve(scene, output, particle=False, polarization='p', one_dimensional=False,
          factorization='default', resolution=8, near_points=(), executable=DEFAULT_EXE,
          timeout=600, epsilon_method='analytic', stripe_rule=False, finite_patch=None):
    output, executable = Path(output), Path(executable).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise FileExistsError('Use an empty output folder to preserve previous raw results')
    if not executable.is_file():
        raise FileNotFoundError('Build the private S4 runtime first: '+str(executable))
    manifest = executable.parent/'build_manifest.json'
    if not manifest.exists():
        manifest = executable.parent.parent/'build_manifest.json'
    build = json.loads(manifest.read_text(encoding='utf-8')) if manifest.exists() else None
    if stripe_rule or scene.harmonics_x is not None:
        if not build or build.get('variant') != 'stripes':
            raise ValueError('Explicit basis/stripe rule requires the separately built stripes runtime')
    source, volume = script(scene, particle, polarization, one_dimensional,
                            factorization, resolution, near_points, epsilon_method, stripe_rule, finite_patch)
    input_path = output/'input.lua'
    input_path.write_text(source, encoding='ascii')
    options = dict(scene=asdict(scene), particle=particle, polarization=polarization,
        one_dimensional=one_dimensional, factorization=factorization, resolution=resolution,
        near_points=list(near_points), epsilon_method=epsilon_method, stripe_rule=stripe_rule)
    if finite_patch is not None:
        options['finite_patch'] = asdict(finite_patch)
    (output/'inputs.json').write_text(json.dumps(options, indent=2), encoding='utf-8')
    env = dict(os.environ, OPENBLAS_NUM_THREADS='2', OMP_NUM_THREADS='2')
    start = time.perf_counter()
    with (output/'stdout.txt').open('w', encoding='utf-8') as stdout, \
            (output/'stderr.txt').open('w', encoding='utf-8') as stderr:
        process = subprocess.run([str(executable), str(input_path.resolve())], env=env,
            stdout=stdout, stderr=stderr, timeout=timeout, check=False)
    if process.returncode:
        raise RuntimeError('S4 failed (0x{:08x}); inspect {}'.format(
            process.returncode & 0xffffffff, output/'stderr.txt'))
    raw = parse((output/'stdout.txt').read_text(encoding='utf-8'))
    result, powers = reconstruct(raw, scene, one_dimensional)
    np.savez_compressed(output/'modes.npz', **result)
    info = dict(options, **powers, actual_modes=len(result['orders']),
        requested_modes=scene.requested_modes,
        retained_order_extent=np.max(abs(result['orders']), axis=0).tolist(),
        absorption_in_finite_stack=1-powers['reflectance']-powers['flux_into_substrate'],
        elapsed_seconds=time.perf_counter()-start, backend='S4', build=build,
        executable_sha256=sha256(executable), input_sha256=sha256(input_path),
        raw_npz_sha256=sha256(output/'modes.npz'),
        geometry=dict(particle_present=particle, sliced_particle_volume_um3=volume,
            exact_particle_volume_um3=np.pi*scene.diameter_um**3/6 if particle else 0),
        coordinates='Global z up; propagation from -KY toward +KY; E phase origin at ridge top z=0.',
        status='EXPERIMENTAL_NOT_CONVERGED',
        scope='Periodic supercell diffraction, not an isolated defect or finite grating.',
        material_caveat='Legacy optical constants are comparison inputs, not a 266 nm calibration.')
    if finite_patch is not None:
        info['geometry']['finite_patch'] = asdict(finite_patch)
        info['scope'] = 'Finite etched motif in a periodic S4 supercell with flat buffers; isolation unverified.'
    (output/'report.json').write_text(json.dumps(info, indent=2), encoding='utf-8')
    return result, info
