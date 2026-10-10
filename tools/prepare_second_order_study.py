"""Prepare bounded second-order haze scans from a previous validated input plan."""
import argparse
import copy
import json
from pathlib import Path

from fullwave_prototype.haze_convergence import validate_plan


def prepare(plan, stage):
    plan = copy.deepcopy(plan)
    p = plan['parameters']
    p.update(order=2, pml_order=1 if stage == 'local' else 2,
             static_condensation=True, direct_solver='sparsecholesky')
    plan['scans'] = (dict(rough_mesh_scale=[1., .9], pml_order=[1, 2]) if stage == 'local'
        else dict(pml_strength=[p['pml_strength'], p['pml_strength']*1.25],
                  pml_thickness_um=[p['pml_thickness_um'], p['pml_thickness_um']*1.125],
                  domain_radius_um=[p['domain_radius_um'], p['domain_radius_um']*1.1]))
    validate_plan(plan)
    return plan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, help='Previous convergence input plan')
    parser.add_argument('--output', required=True)
    parser.add_argument('--stage', choices=('local', 'boundary'), required=True)
    args = parser.parse_args()
    plan = prepare(json.loads(Path(args.input).read_text(encoding='utf-8')), args.stage)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(plan, indent=2), encoding='utf-8')
    print(str(output))


if __name__ == '__main__':
    main()
