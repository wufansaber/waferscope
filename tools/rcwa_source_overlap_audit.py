"""Independent real-space polar quadrature check of stored particle Fourier fields."""
import argparse
import json
from pathlib import Path

import numpy as np

from rcwa_prototype.current_overlap import compare_source_fields, validate_source
from rcwa_prototype.current_radiation import load_source


def polar_overlap(first, second, radial, azimuthal):
    a, ka, za, weights, radii = validate_source(first)
    b, kb, zb, wb, rb = validate_source(second)
    for x, y in ((za, zb), (weights, wb), (radii, rb)):
        np.testing.assert_allclose(x, y, rtol=2e-13, atol=1e-15)
    nodes, wr = np.polynomial.legendre.leggauss(radial)
    phi = np.arange(azimuthal)*2*np.pi/azimuthal
    norm, difference = 0., 0.
    for radius in np.unique(radii):
        take = radii == radius
        r = (nodes+1)*radius/2
        rr, pp = np.meshgrid(r, phi, indexing='ij')
        points = np.column_stack((rr.ravel()*np.cos(pp.ravel()), rr.ravel()*np.sin(pp.ravel())))
        w = np.repeat(wr*radius/2*r, azimuthal)*2*np.pi/azimuthal

        def field(coefficients, wavevectors):
            value = np.exp(1j*points@wavevectors.T)@coefficients[take].transpose(1, 0, 2).reshape(len(wavevectors), -1)
            return value.reshape(len(points), int(take.sum()), 3)

        ea, eb = field(a, ka), field(b, kb)
        norm += np.einsum('pz,p,z->', np.sum(abs(ea)**2, axis=-1), w, weights[take])
        difference += np.einsum('pz,p,z->', np.sum(abs(eb-ea)**2, axis=-1), w, weights[take])
    return dict(radial=radial, azimuthal=azimuthal,
        first_internal_e_volume_norm_um3=float(norm), difference_volume_norm_um3=float(difference),
        internal_e_relative_l2=float(np.sqrt(difference/norm)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first', type=Path, required=True)
    parser.add_argument('--second', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Use a new output file to preserve prior audit results')
    a, scene, _, pa = load_source(args.first)
    b, _, _, pb = load_source(args.second)
    exact = compare_source_fields(a, b, scene.eps_particle-1)
    levels = [polar_overlap(a, b, n, m) for n, m in ((20, 80), (32, 128))]
    for row in levels:
        row['relative_l2_change_vs_analytic_disk'] = row['internal_e_relative_l2']/exact['internal_e_relative_l2']-1
    result = dict(first=pa, second=pb, analytic_disk=exact, polar_quadrature=levels,
        scope='Same stored truncated fields. This checks the overlap integration, not Maxwell or isolated-defect accuracy.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
