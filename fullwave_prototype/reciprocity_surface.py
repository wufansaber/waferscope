"""Independent Lorentz-reciprocity surface extraction above a planar half-space.

Uses scattered E and H*Z0 on a closed box outside the entire material contrast.
No complex conjugation in the reciprocal surface flux (not a power integral).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss


def box_quadrature(bounds, order=48):
    bounds = np.asarray(bounds, float)
    if bounds.shape != (3, 2) or not np.isfinite(bounds).all() or np.any(bounds[:, 1] <= bounds[:, 0]):
        raise ValueError('Three increasing finite box intervals required')
    if order < 4:
        raise ValueError('At least four Gauss nodes per interval required')
    nodes, weights = leggauss(order)
    points, normals, areas = [], [], []
    for fixed in range(3):
        moving = [d for d in range(3) if d != fixed]
        intervals = []
        for d in moving:
            a, b = bounds[d]
            intervals.append([(a, 0), (0, b)] if d == 2 and a < 0 < b else [(a, b)])
        for side in (0, 1):
            for lo0, hi0 in intervals[0]:
                for lo1, hi1 in intervals[1]:
                    v0, v1 = (nodes+1)*(hi0-lo0)/2+lo0, (nodes+1)*(hi1-lo1)/2+lo1
                    xx, yy = np.meshgrid(v0, v1)
                    block = np.zeros((xx.size, 3))
                    block[:, fixed] = bounds[fixed, side]
                    block[:, moving[0]], block[:, moving[1]] = xx.ravel(), yy.ravel()
                    normal = np.zeros_like(block)
                    normal[:, fixed] = -1 if side == 0 else 1
                    points.append(block)
                    normals.append(normal)
                    areas.append((np.outer(weights, weights)*(hi0-lo0)*(hi1-lo1)/4).ravel())
    return np.concatenate(points), np.concatenate(normals), np.concatenate(areas)


def amplitude(points, normals, areas, electric, magnetic, directions, p):
    from .planar import fresnel
    points, directions = np.asarray(points), np.asarray(directions)
    if (points.shape != np.shape(electric) or points.shape != np.shape(magnetic) or
            points.shape != np.shape(normals) or np.shape(areas) != (len(points),) or
            directions.ndim != 2 or directions.shape[1] != 3):
        raise ValueError('Mismatching surface fields, normals or directions')
    if not all(np.isfinite(v).all() for v in (points, normals, areas, electric, magnetic, directions)):
        raise ValueError('Finite reciprocity inputs required')
    if np.any(directions[:, 2] <= 0) or not np.allclose(np.sum(directions**2, axis=1), 1):
        raise ValueError('Unit upper-hemisphere observation directions required')
    output = np.zeros((len(directions), 2), complex)
    above = points[:, 2] >= 0
    weighted_normal = np.asarray(normals)*np.asarray(areas)[:, None]
    for i, k in enumerate(directions):
        rho = np.hypot(k[0], k[1])
        s = np.array((-k[1]/rho, k[0]/rho, 0)) if rho > 1e-12 else np.array((-1., 0, 0))
        down = k*np.array((1, 1, -1))
        rs, rp, q = fresnel(k[2], p.eps_substrate)
        for pol, incident in enumerate((s, np.cross(k, s))):
            reflected = rs*s if pol == 0 else rp*np.cross(down, s)
            transmitted = np.array(incident+reflected, complex)
            transmitted[2] = -(k[0]*transmitted[0]+k[1]*transmitted[1])/q
            kt = np.array((-k[0], -k[1], -q))
            direct_phase = np.exp(-1j*p.k0*(points@k))
            reflected_phase = np.exp(-1j*p.k0*(points@down))
            sub_phase = np.exp(1j*p.k0*(points@kt))
            et = np.where(above[:, None], direct_phase[:, None]*incident+
                reflected_phase[:, None]*reflected, sub_phase[:, None]*transmitted)
            ht = np.where(above[:, None], direct_phase[:, None]*np.cross(-k, incident)+
                reflected_phase[:, None]*np.cross(-down, reflected),
                sub_phase[:, None]*np.cross(kt, transmitted))
            flux = np.cross(electric, ht)-np.cross(et, magnetic)
            output[i, pol] = 1j*p.k0/(4*np.pi)*np.sum(flux*weighted_normal)
    return output


def evaluate_half_field(mesh, field, points, p):
    import ngsolve as ng
    mapped = np.array(points, copy=True)
    negative = mapped[:, 0] < 0
    mapped[:, 0] = abs(mapped[:, 0])
    location = mesh(mapped[:, 0], mapped[:, 1], mapped[:, 2])
    e = np.asarray(field(location)).reshape(-1, 3)
    h = np.asarray(ng.curl(field)(location)).reshape(-1, 3)/(1j*p.k0)
    e[negative] *= (-1, 1, 1)
    h[negative] *= (1, -1, -1)
    return e, h
