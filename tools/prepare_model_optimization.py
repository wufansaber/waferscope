"""Reproducible resource-bounded order, coupled-height and angle validation plan."""
import argparse
import json
from pathlib import Path

from fullwave_prototype.haze_convergence import default_plan, validate_plan
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.rough_surface import Roughness


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    p = Parameters(wavelength_nm=266, gap_nm=10, order=1, threads=2,
                   direct_solver='sparsecholesky', static_condensation=True,
                   pml_order=1, max_dofs=500000)
    r = Roughness(model='abc', grid_points=9, taper_fraction=.125,
                  correlation_length_um=.2, psd_exponent=4.5)
    plan = default_plan(p, r, samples=3)
    plan['scans'] = dict(fem_order=[1, 2], height_grid_points=[9, 13, 17],
                         angular_step_deg=[10, 5, 2.5])
    validate_plan(plan)
    (output/'optimization_plan.json').write_text(json.dumps(plan, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
