"""Bounded-memory electric-only sampling of the existing coherent Fresnel rays."""
import numpy as np

from .illumination import illumination_rays


def electric_background(points, p, block=512):
    return sampled_background(points, p, block)[0]


def sampled_background(points, p, block=512, magnetic=False):
    points = np.asarray(points, float)
    if (points.ndim < 1 or points.shape[-1] != 3 or not np.isfinite(points).all() or
            not isinstance(block, int) or block < 1):
        raise ValueError('Finite XYZ points and positive sampling block required')
    flat = points.reshape(-1, 3)
    rays = illumination_rays(p)
    down, up, kt, ei, er, et = (np.array([ray[i] for ray in rays]) for i in range(6))
    result = np.empty(flat.shape, complex)
    hresult = np.empty(flat.shape, complex) if magnetic else None
    if magnetic:
        hi, hr, ht = np.cross(down, ei), np.cross(up, er), np.cross(kt, et)
    for start in range(0, len(flat), block):
        xyz = flat[start:start+block]
        air = xyz[:, 2] >= 0
        values = np.empty(xyz.shape, complex)
        if air.any():
            pi = np.exp(1j*p.k0*(xyz[air]@down.T))
            pr = np.exp(1j*p.k0*(xyz[air]@up.T))
            values[air] = pi@ei+pr@er
            if magnetic:
                hresult[start:start+len(xyz)][air] = pi@hi+pr@hr
        if (~air).any():
            pt = np.exp(1j*p.k0*(xyz[~air]@kt.T))
            values[~air] = pt@et
            if magnetic:
                hresult[start:start+len(xyz)][~air] = pt@ht
        result[start:start+block] = values
    return result.reshape(points.shape), hresult.reshape(points.shape) if magnetic else None
