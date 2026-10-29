"""Small native-candidate checks, not an AFMM defect or imaging certification."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

from tools.build_crcwa_runtime import COMMIT


def check_interface(data):
    # Upstream `for :indexyz 1 10 1` has an exclusive upper bound.
    indices = np.arange(1, 10)
    if data.shape != (9, 3) or not np.isfinite(data).all():
        raise ValueError('Expected all nine upstream interface angles')
    theta = np.arcsin(indices*700e-9/(1.5*5e-6))
    cost = np.sqrt(1-(1.5/5.)**2*np.sin(theta)**2)
    expected = abs((1.5*np.cos(theta)-5*cost)/(1.5*np.cos(theta)+5*cost))**2
    absolute_error = float(np.max(abs(data[:, 1]-expected)))
    angle_error = float(np.max(abs(data[:, 0]-np.rad2deg(theta))))
    if absolute_error > 2e-5 or angle_error > 2e-4:
        raise RuntimeError('Independent Fresnel or angle check did not pass')
    return theta, expected, absolute_error, angle_error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', type=Path, default=Path('build/crcwa_runtime_20261027'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    runtime, output = args.runtime.resolve(), args.output.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error('Use an empty directory to retain candidate benchmark records')
    source = runtime/('cRCWA-'+COMMIT)
    executable = runtime/'bin/crcwa.exe'
    manifest = json.loads((runtime/'build_manifest.json').read_text())
    if hashlib.sha256(executable.read_bytes()).hexdigest() != manifest['executable_sha256']:
        parser.error('Candidate executable hash does not match its build manifest')
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    output.mkdir(parents=True, exist_ok=True)
    records = []
    for name, script in (('plane_interface', source/'test/1D_test_angle/1D_test_angle.fmm'),
                         ('upstream_simple_mode', source/'test/simplemode/simplemode.fmm')):
        folder = output/name
        folder.mkdir()
        shutil.copy2(script, folder/script.name)
        with (folder/'stdout.txt').open('w', encoding='utf-8') as stream:
            result = subprocess.run([str(executable), script.name], cwd=folder, env=env,
                stdout=stream, stderr=subprocess.STDOUT, timeout=90)
        log = (folder/'stdout.txt').read_text(encoding='utf-8', errors='replace')
        if (result.returncode or not (folder/'test_results.txt').is_file()
                or 'Could not recognize the command' in log):
            raise RuntimeError('Native candidate check failed: '+name)
        records.append(dict(name=name, script_sha256=hashlib.sha256(script.read_bytes()).hexdigest(),
                            returncode=result.returncode))
    data = np.loadtxt(output/'plane_interface/test_results.txt')
    theta, expected, absolute_error, angle_error = check_interface(data)
    original = (source/'test/simplemode/test_reference.txt').read_text().split()
    actual = (output/'upstream_simple_mode/test_results.txt').read_text().split()
    reference_neff, observed_neff = float(original[2]), float(actual[2])
    if abs(reference_neff-observed_neff) > 1e-5:
        raise RuntimeError('Upstream printed mode regression did not pass')
    np.savez_compressed(output/'plane_interface_raw.npz', angle_deg=data[:, 0],
        calculated_reflectance=data[:, 1], upstream_script_fresnel=data[:, 2],
        independent_fresnel=expected, independent_angle_deg=np.rad2deg(theta))
    report = dict(source_commit=COMMIT, build_manifest=str(runtime/'build_manifest.json'),
        executable_sha256=manifest['executable_sha256'], runs=records,
        plane_interface=dict(wavelength_nm=700., incident_n=1.5, transmitted_n=5.,
            polarization='s (Ex with incidence in YZ)', angles=int(len(theta)),
            angle_loop='Upstream exclusive range: indices 1 through 9.',
            max_absolute_reflectance_error=absolute_error, max_angle_error_deg=angle_error,
            output_precision_note='Upstream text output has about six significant digits; not a full-double error bound.'),
        upstream_mode=dict(printed_reference_neff=reference_neff, printed_observed_neff=observed_neff,
            scope='Upstream small periodic waveguide regression, not independent AFMM/PML validation.'),
        status='SMALL_INTERFACE_AND_UPSTREAM_REGRESSION_PASSED_NOT_AFMM_ACCEPTANCE',
        limitations=['No PML isolation, high-index grating plus sphere, near-field continuity or imaging acceptance yet.',
                     'Two Windows compatibility patches: MinGW semaphore header and binary script offsets; no Maxwell formula changes.',
                     'Candidate remains outside the application package and GUI.'])
    (output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
