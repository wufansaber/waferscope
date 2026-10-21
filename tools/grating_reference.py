"""NIST RCW oracle; run under NIST2, separately from native NGSolve."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from fullwave_prototype.periodic_grating import GratingInputs


def main():
    from pySCATMECH.rcw import RCW_Model
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--order', type=int, default=45)
    parser.add_argument('--width', type=float, default=.25)
    args = parser.parse_args()
    if args.order < 12:
        parser.error('At least 12 Fourier orders are required for the exported index range')
    p = GratingInputs(ridge_width_um=args.width)
    p.validate()
    grating = {None: 'Single_Line_Grating', 'period': p.period_um,
        'medium_i': 1, 'medium_t': complex(p.substrate_n, p.substrate_k),
        'material': complex(p.substrate_n, p.substrate_k), 'space': 1,
        'height': p.depth_um, 'topwidth': p.ridge_width_um,
        'bottomwidth': p.ridge_width_um, 'offset': 0, 'nlevels': 1}
    model = RCW_Model(dict(order=args.order, type=0, **{'lambda': p.wavelength_um,
        'thetai': p.incidence_deg, 'rotation': 0, 'grating': grating}))
    indices = np.arange(-12, 13)
    direction, amplitude, efficiency = [], [], []
    for m in indices:
        jones = np.asarray(model.DiffractionAmplitude(int(-m)))
        mueller = np.asarray(model.DiffractionEfficiency(int(-m)))
        direction.append(model.Direction(int(-m)))
        amplitude.append(jones[1, 1])
        efficiency.append(mueller[0, 0]-mueller[0, 1])
    # Record actual NIST geometry instead of assuming its lateral fiducial.
    x = np.linspace(-.249, .249, 101)
    eps = np.array([model.getEpsilon(float(v), -p.depth_um/2) for v in x])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, orders=indices, direction=direction,
        amplitude=amplitude, efficiency=efficiency, geometry_x_um=x, geometry_epsilon=eps)
    args.output.with_suffix('.json').write_text(json.dumps(dict(inputs=asdict(p),
        fourier_order=args.order, polarization='p', particle=False,
        order_mapping='canonical m=-NIST index; canonical transverse k/k0=sin(theta)+m*lambda/period',
        reflectance=float(np.sum(efficiency)),
        scope='Infinite 1D rectangular silicon ridges, invariant along the lines; not a finite SEM reconstruction.'),
        indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
