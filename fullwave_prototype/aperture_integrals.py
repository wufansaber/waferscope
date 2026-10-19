"""Solid-angle quadrature of circular objectives, clipped at the wafer horizon."""
from functools import lru_cache

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RegularGridInterpolator


@lru_cache(maxsize=32)
def cap_quadrature(na, tilt_degrees, order=20):
    if (not np.isfinite(na) or not 0 <= na < 1 or not np.isfinite(tilt_degrees) or
            abs(tilt_degrees) > 89.9 or type(order) is not int or not 2 <= order <= 256):
        raise ValueError('Cap requires NA [0,1), tilt [-89.9,89.9], quadrature order [2,256]')
    if na == 0:
        return np.empty((0, 3)), np.empty(0)
    tilt = np.deg2rad(tilt_degrees)
    sine, cosine = np.sin(tilt), np.cos(tilt)
    lower_cap = np.sqrt(1-na*na)
    breaks = [0., 2*np.pi]
    if abs(sine) > 1e-15:
        threshold = cosine*lower_cap/(sine*na)
        if abs(threshold) < 1:
            angle = np.arccos(threshold)
            breaks += [angle, 2*np.pi-angle]
    nodes, weights = leggauss(order)
    directions, solid_angles = [], []
    axis, transverse = np.array((sine, 0, cosine)), np.array((cosine, 0, -sine))
    for start, end in zip(sorted(breaks)[:-1], sorted(breaks)[1:]):
        beta = (nodes+1)*(end-start)/2+start
        beta_weight = weights*(end-start)/2
        a = np.maximum(sine*np.cos(beta), 0.)
        horizon = a/np.hypot(a, cosine)
        lower = np.maximum(lower_cap, horizon)
        span = np.where(horizon > lower_cap, 1-horizon, na*na/(1+lower_cap))
        mu = lower[:, None]+(nodes[None, :]+1)*span[:, None]/2
        radial = np.sqrt(np.maximum(0, (1-mu)*(1+mu)))
        ray = (mu[..., None]*axis + radial[..., None]*(
            np.cos(beta)[:, None, None]*transverse +
            np.sin(beta)[:, None, None]*np.array((0., 1., 0.))))
        directions.append(ray.reshape(-1, 3))
        solid_angles.append((beta_weight[:, None]*weights[None, :]*span[:, None]/2).ravel())
    directions, solid_angles = np.concatenate(directions), np.concatenate(solid_angles)
    directions.setflags(write=False)
    solid_angles.setflags(write=False)
    return directions, solid_angles


def integrate_maps(maps, theta, phi, apertures, tilts, order=20):
    """Integrate periodic, bilinearly interpolated DSC, not a cap's geometric area."""
    values = np.moveaxis(np.asarray(maps), 0, -1)
    periodic_values = np.concatenate((values, values[:, :1]), axis=1)
    periodic_phi = np.concatenate((phi, [phi[0]+2*np.pi]))
    interpolate = RegularGridInterpolator((theta, periodic_phi), periodic_values,
                                         bounds_error=True)
    totals = np.zeros((values.shape[-1], len(apertures)))
    for index, (na, tilt) in enumerate(zip(apertures, tilts)):
        directions, weights = cap_quadrature(float(na), float(tilt), order)
        if not len(weights):
            continue
        polar = np.arccos(np.clip(directions[:, 2], 0, 1))
        azimuth = (np.arctan2(directions[:, 1], directions[:, 0])-phi[0]) % (2*np.pi)+phi[0]
        samples = interpolate(np.column_stack((polar, azimuth)))
        totals[:, index] = weights @ samples
    return totals


def integrate_sources(points, sources, parameters, apertures, tilts, order=20):
    """Independent reference: evaluate the saved complex FEM volume sources directly."""
    from .far_field import radiation_amplitude
    totals = []
    for na, tilt in zip(apertures, tilts):
        directions, weights = cap_quadrature(float(na), float(tilt), order)
        if not len(weights):
            totals.append(0.)
            continue
        amplitude = radiation_amplitude(points, sources, directions, parameters)
        totals.append(float(weights @ np.sum(abs(amplitude)**2, axis=1)))
    return np.array(totals)
