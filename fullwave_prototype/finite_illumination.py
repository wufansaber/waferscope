"""Propagating coherent angular spectrum of an elliptical wafer Gaussian."""
from functools import lru_cache

import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss


def gaussian_spectrum(p, order=None):
    """Amplitude quadrature: X/Y target intensity 1/e^2 radii at z=0."""
    order = order or p.beam_order
    if order > 151:
        from scipy.special import roots_hermite, roots_legendre
        nodes, weights = roots_hermite(order)
        gy, wy = roots_legendre(order)
    else:
        nodes, weights = hermgauss(order)
        gy, wy = leggauss(order)
    wx, length = p.beam_width_um/2, p.beam_length_um/2
    center = np.sin(np.deg2rad(p.incidence_deg))
    directions, amplitudes = [], []
    for x, weight in zip(nodes, weights):
        qx = 2*x/(p.k0*wx)
        if abs(qx) >= 1:
            continue
        edge = np.sqrt(1-qx*qx)
        low = max(-6., (-edge-center)*p.k0*length/2)
        high = min(6., (edge-center)*p.k0*length/2)
        if high <= low:
            continue
        if low == -6. and high == 6.:
            y, y_weights = nodes, weights
        else:
            # Smooth the grazing-wave endpoint without undersampling a full Gaussian.
            y = (high+low)/2+(high-low)/2*np.sin(np.pi*gy/2)
            jacobian = (high-low)*np.pi/4*np.cos(np.pi*gy/2)
            y_weights = wy*jacobian*np.exp(-y*y)
        qy = center+2*y/(p.k0*length)
        qz = -np.sqrt(np.maximum(0, 1-qx*qx-qy*qy))
        a = weight*y_weights/np.pi
        keep = a > 1e-20 if order > 151 else np.ones(order, bool)
        if np.any(keep):
            directions.append(np.column_stack((np.full(order, qx), qy, qz))[keep])
            amplitudes.append(a[keep])
    if not directions:
        raise ValueError('Gaussian angular spectrum has no propagating quadrature rays')
    return np.concatenate(directions), np.concatenate(amplitudes)


def _raw_power_area(p, order, flux=True):
    """Parseval integral of downward Poynting flux, not rectangular P/A."""
    nodes, weights = leggauss(order)
    wx, wy = p.beam_width_um/2, p.beam_length_um/2
    center = np.sin(np.deg2rad(p.incidence_deg))
    bound = min(6., p.k0*wx/2)
    x = nodes*bound
    qx = 2*x/(p.k0*wx)
    edge = np.sqrt(np.maximum(0, 1-qx*qx))
    low = np.maximum(-6., (-edge-center)*p.k0*wy/2)
    high = np.minimum(6., (edge-center)*p.k0*wy/2)
    span = np.maximum(0, high-low)
    y = (high+low)[:, None]/2+span[:, None]/2*np.sin(np.pi*nodes/2)
    jacobian = span[:, None]*np.pi/4*np.cos(np.pi*nodes/2)
    qy = center+2*y/(p.k0*wy)
    cosine = np.sqrt(np.maximum(0, 1-qx[:, None]**2-qy*qy))
    integrand = np.exp(-2*(x[:, None]**2+y*y))*(cosine if flux else 1.)
    return float(wx*wy*bound*np.sum(weights[:, None]*weights*jacobian*integrand))


@lru_cache(maxsize=64)
def power_area(p):
    from .illumination import raw_rays, central_flux
    from .beam_inputs import raw_area, wafer_dimensions
    from dataclasses import replace
    length, width = wafer_dimensions(p)
    effective = replace(p, beam_length_um=length, beam_width_um=width)
    integrate = (lambda order: _raw_power_area(effective, order)) if p.beam_definition in (
        'wafer_gaussian', 'beam_cross_section') else (lambda order: raw_area(p, order))
    value = integrate(192)/central_flux(p, raw_rays(p))
    check = integrate(128)/central_flux(p, raw_rays(p))
    if not np.isfinite(value) or value <= 0 or abs(check/value-1) > 1e-5:
        raise ValueError('Finite beam power quadrature did not converge; revise beam dimensions')
    return value


def incident_fields(points, p):
    from .illumination import illumination_rays
    points = np.asarray(points)
    e, h = np.zeros(points.shape, complex), np.zeros(points.shape, complex)
    for down, _, _, ei, _, _ in illumination_rays(p):
        phase = np.exp(1j*p.k0*(points@down))[..., None]
        e += phase*ei
        h += phase*np.cross(down, ei)
    return e, h


def patch_coverage(p, width_um, order=40):
    nodes, weights = leggauss(order)
    x, y = np.meshgrid(nodes*width_um/2, nodes*width_um/2, indexing='xy')
    points = np.stack((x, y, np.zeros_like(x)), axis=-1)
    e, h = incident_fields(points, p)
    flux = -np.real(np.cross(e, h.conj()))[..., 2]
    integral = float((width_um/2)**2*np.sum(weights[:, None]*weights*flux))
    return integral/power_area(p)
