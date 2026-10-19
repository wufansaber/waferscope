"""Memory-bounded matrix evaluation of half-space volume radiation."""
import numpy as np

from .planar import fresnel


def radiation_amplitude(points, weights, directions, p, block_size=32):
    points, weights, directions = np.asarray(points), np.asarray(weights), np.asarray(directions)
    if (points.ndim != 2 or points.shape[1] != 3 or weights.shape != points.shape
            or directions.ndim != 2 or directions.shape[1] != 3
            or not isinstance(block_size, int) or block_size < 1):
        raise ValueError('Expected N by 3 volume sources and directions, positive block size')
    if not all(np.isfinite(a).all() for a in (points, weights, directions)):
        raise ValueError('Non-finite radiation input')
    if ((p.eps_substrate != 1 and np.any(directions[:, 2] < -1e-12))
            or not np.allclose(np.sum(directions**2, axis=1), 1)):
        raise ValueError('Unit directions required; a material substrate supports upper hemisphere only')
    result = np.zeros((len(directions), 2), complex)
    above = points[:, 2] >= 0
    up, down, wu, wd = points[above], points[~above], weights[above], weights[~above]
    # Limit direction blocking; one direction still requires N_source entries.
    block_size = min(block_size, max(1, 2000000//max(1, len(points))))
    for start in range(0, len(directions), block_size):
        k = directions[start:start+block_size]
        rho = np.hypot(k[:, 0], k[:, 1])
        s = np.column_stack((-k[:, 1], k[:, 0], np.zeros(len(k))))/np.maximum(rho[:, None], 1e-15)
        s[rho < 1e-12] = (-1, 0, 0)
        kd = k*np.array((1, 1, -1))
        pup, pd = np.cross(k, s), np.cross(kd, s)
        rs, rp, q = fresnel(k[:, 2], p.eps_substrate)
        if p.eps_substrate == 1:
            q = k[:, 2]
        direct = np.exp(-1j*p.k0*(k@up.T))@wu
        reflected = np.exp(-1j*p.k0*(kd@up.T))@wu
        value = np.column_stack((np.sum(s*(direct+rs[:, None]*reflected), axis=1),
            np.sum(pup*direct+rp[:, None]*pd*reflected, axis=1)))
        if len(down):
            ts = (1+rs[:, None])*s
            tp = np.asarray((1-rp[:, None])*pup, complex)
            if p.eps_substrate != 1:
                tp[:, 2] = -2*np.sum(k[:, :2]*pup[:, :2], axis=1)/(p.eps_substrate*k[:, 2]+q)
            kt = np.column_stack((k[:, :2], q))
            source = np.exp(-1j*p.k0*(kt@down.T))@wd
            value += np.column_stack((np.sum(ts*source, axis=1), np.sum(tp*source, axis=1)))
        result[start:start+len(k)] = p.k0**2/(4*np.pi)*value
    return result
