"""Homogeneous-air surface-current far field, independent of volume contrast sources."""
import argparse
import json
from pathlib import Path

import numpy as np

from .parameters import Parameters
from .planar import background_fields


def surface_amplitude(points, normals, weights, electric, hz0, directions, k0):
    js = np.cross(normals, hz0)*weights[:, None]
    ms = -np.cross(normals, electric)*weights[:, None]
    result = []
    for start in range(0, len(directions), 64):
        k = directions[start:start+64]
        phase = np.exp(-1j*k0*(k@points.T))
        j, m = phase@js, phase@ms
        result.append(1j*k0/(4*np.pi)*(j-k*np.sum(k*j, axis=1)[:, None]-np.cross(k, m)))
    return np.concatenate(result)


def check_saved(directory):
    from .physical_audit import sphere_rule
    directory = Path(directory)
    report = json.loads((directory/'report.json').read_text(encoding='utf-8'))
    p = Parameters(**report['parameters'])
    if p.eps_substrate != 1:
        raise ValueError('This Huygens test requires a homogeneous air exterior')
    with np.load(directory/'fields.npz') as data:
        arrays = dict(data)
    directions, measure = sphere_rule(20)
    np.testing.assert_allclose(directions, arrays['directions'], atol=0, rtol=0)
    rho = np.hypot(directions[:, 0], directions[:, 1])
    s = np.column_stack((-directions[:, 1]/rho, directions[:, 0]/rho, np.zeros(len(rho))))
    cart_volume = (s*arrays['amplitude'][:, 0, None]+
                   np.cross(directions, s)*arrays['amplitude'][:, 1, None])
    comparisons, exact_shells = [], []
    count = len(directions)
    for i in (2, 3):
        radius = report['near_shells'][i]['radius_nm']/1000
        sl = slice(i*count, (i+1)*count)
        points = arrays['points_um'][sl]
        eb, hb = background_fields(points, p)
        fem = surface_amplitude(points, directions, measure*radius**2,
            arrays['e'][sl]-eb, arrays['hz0'][sl]-hb, directions, p.k0)
        exact = surface_amplitude(points, directions, measure*radius**2,
            arrays['exact_e'][sl]-eb, arrays['exact_hz0'][sl]-hb, directions, p.k0)
        weight = np.sqrt(measure)[:, None]
        norm = np.linalg.norm(exact*weight)
        comparisons.append(dict(radius_nm=radius*1000,
            surface_vs_volume_complex_l2=float(np.linalg.norm((fem-cart_volume)*weight)/norm),
            surface_vs_mie_complex_l2=float(np.linalg.norm((fem-exact)*weight)/norm),
            volume_vs_mie_complex_l2=float(np.linalg.norm((cart_volume-exact)*weight)/norm),
            exact_surface_csca_relative_error=abs(float(np.sum(measure*np.sum(abs(exact)**2, axis=1)))/report['exact']['csca_um2']-1)))
        exact_shells.append(exact)
    output = dict(comparisons=comparisons,
        analytic_surface_radius_invariance=float(np.linalg.norm(exact_shells[0]-exact_shells[1])/np.linalg.norm(exact_shells[1])),
        note='Air only; enclosing-surface E/H equivalence independent of volume extractor. FEM surface angular quadrature remains a separate discretization error.')
    (directory/'huygens.json').write_text(json.dumps(output, indent=2, allow_nan=False), encoding='utf-8')
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directories', nargs='+')
    args = parser.parse_args()
    for directory in args.directories:
        print(directory, json.dumps(check_saved(directory), indent=2), flush=True)


if __name__ == '__main__':
    main()
