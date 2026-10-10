"""Finite angular-spectrum illumination, with exact Fresnel matching per ray."""
import numpy as np
from functools import lru_cache
from numpy.polynomial.hermite import hermgauss

from .planar import background_vectors, fresnel


@lru_cache(maxsize=64)
def raw_rays(p):
    if p.beam_profile == 'gaussian':
        from .beam_inputs import spectrum
        down0, _, _, e0, _, _ = background_vectors(p)
        directions, weights = spectrum(p)
        return [matched_ray(p, down0, e0, down, weight)
                for down, weight in zip(directions, weights)]
    if p.beam_waist_um == 0:
        return [background_vectors(p)]
    down0, _, _, e0, _, _ = background_vectors(p)
    s0 = np.array((-1., 0, 0))
    p0 = np.cross(down0, s0)
    nodes, weights = hermgauss(p.beam_order)
    rays = []
    for i, x in enumerate(nodes):
        for j, y in enumerate(nodes):
            u, v = 2*np.array((x, y))/(p.k0*p.beam_waist_um)
            if u*u+v*v >= 1:
                raise ValueError('Beam angular quadrature includes evanescent waves; increase waist.')
            down = u*s0+v*p0+np.sqrt(1-u*u-v*v)*down0
            if down[2] >= -1e-10:
                raise ValueError('Beam angular quadrature crosses the wafer horizon; increase waist.')
            axis = np.cross(down0, down)
            ei = e0+np.cross(axis, e0)+np.cross(axis, np.cross(axis, e0))/(1+down0@down)
            rho = np.linalg.norm(down[:2])
            s = np.array((-down[1], down[0], 0))/rho if rho > 1e-14 else s0
            pd = np.cross(down, s)
            up = down*np.array((1, 1, -1))
            c = -down[2]
            rs, rp, q = fresnel(c, p.eps_substrate)
            kt = np.array((down[0], down[1], -q), complex)
            a, b = ei@s, ei@pd
            er = a*rs*s+b*rp*np.cross(up, s)
            tp = (1-rp)*pd.astype(complex)
            tp[2] = 2*rho*c/(p.eps_substrate*c+q)
            et = a*(1+rs)*s+b*tp
            weight = weights[i]*weights[j]/np.pi
            rays.append((down, up, kt, weight*ei, weight*er, weight*et))
    return rays


def matched_ray(p, down0, e0, down, weight):
    axis = np.cross(down0, down)
    ei = e0+np.cross(axis, e0)+np.cross(axis, np.cross(axis, e0))/(1+down0@down)
    rho = np.linalg.norm(down[:2])
    s = np.array((-down[1], down[0], 0))/rho if rho > 1e-14 else np.array((-1., 0, 0))
    pd, up = np.cross(down, s), down*np.array((1, 1, -1))
    c = -down[2]
    rs, rp, q = fresnel(c, p.eps_substrate)
    kt = np.array((down[0], down[1], -q), complex)
    a, b = ei@s, ei@pd
    er = a*rs*s+b*rp*np.cross(up, s)
    tp = (1-rp)*pd.astype(complex)
    tp[2] = 2*rho*c/(p.eps_substrate*c+q)
    et = a*(1+rs)*s+b*tp
    return down, up, kt, weight*ei, weight*er, weight*et


def central_flux(p, rays):
    if p.beam_profile == 'gaussian' and p.beam_definition in ('wafer_rectangle', 'measured_wafer'):
        # Unit input-envelope peak E, not central irradiance (which may vanish).
        return 1.
    down0 = background_vectors(p)[0]
    e = sum(ray[3] for ray in rays)
    h = sum(np.cross(ray[0], ray[3]) for ray in rays)
    value = float(np.real(np.cross(e, h.conj()))@down0)
    if not np.isfinite(value) or value <= 0:
        raise ValueError('Incident central flux must be positive and finite')
    return value


@lru_cache(maxsize=64)
def illumination_rays(p):
    rays = raw_rays(p)
    scale = np.sqrt(central_flux(p, rays))
    return [tuple(ray[:3])+tuple(v/scale for v in ray[3:]) for ray in rays]
