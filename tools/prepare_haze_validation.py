"""Synthetic validation fixtures, explicitly not measured wafer data."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from fullwave_prototype.experimental_input import prepare
from fullwave_prototype.haze_convergence import default_plan, validate_plan
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.rough_surface import Roughness


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    p = Parameters(wavelength_nm=266, gap_nm=10, order=1)
    r = Roughness(model='abc', correlation_length_um=.2, psd_exponent=4.5,
                  grid_points=9, taper_fraction=.125)
    plan = default_plan(p, r, samples=4)
    plan['scans'].update(patch_width_um=[.2, .225, .25],
                        taper_width_um=[.025, .05, .075])
    validate_plan(plan)
    (output/'convergence_plan.json').write_text(json.dumps(plan, indent=2), encoding='utf-8')
    sparse_plan = dict(plan, parameters=dict(asdict(p), direct_solver='sparsecholesky'))
    (output/'convergence_sparse_plan.json').write_text(json.dumps(sparse_plan, indent=2), encoding='utf-8')
    sparse_parameters = dict(asdict(p), direct_solver='sparsecholesky')
    sparse_roughness = dict(asdict(r), seed=9, grid_points=13)
    (output/'sparse_solver_check.json').write_text(json.dumps(dict(
        parameters=sparse_parameters, roughness=sparse_roughness,
        angular_step_deg=5, quadrature_order=5), indent=2), encoding='utf-8')
    axis = np.arange(17)*.0125
    x, y = np.meshgrid(axis, axis)
    height = 1.8*np.cos(2*np.pi*x/.1)*np.cos(2*np.pi*y/.08)+3*x-2*y
    source = output/'synthetic_height.npz'
    np.savez_compressed(source, height_nm=height, dx_um=.0125, dy_um=.0125)
    source.with_suffix('.json').write_text(json.dumps(dict(
        origin='Synthetic sinusoidal validation fixture, not experimental wafer data',
        instrument='analytic generator', calibration='nm and um specified explicitly')), encoding='utf-8')
    report = prepare(source, output/'measurement_audit')
    for model, path in (('measured_height', source),
                        ('measured_psd', output/'measurement_audit'/'measured_psd.npz')):
        measured = dict(asdict(r), model=model, measurement_path=str(path), rms_nm=report['rms_nm'])
        measured['grid_points'] = 9
        inputs = dict(parameters=asdict(p), roughness=measured, angular_step_deg=5, quadrature_order=5)
        (output/(model+'_input.json')).write_text(json.dumps(inputs, indent=2), encoding='utf-8')
    print(output)


if __name__ == '__main__':
    main()
