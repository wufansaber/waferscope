"""Integrate the complex surface radiation directly, independently of its map."""
import numpy as np
from numpy.polynomial.legendre import leggauss

from .aperture_integrals import cap_quadrature


def integrals(radiation, apertures, tilts, order):
    nodes, weights = leggauss(order)
    mu, weight_mu = (nodes+1)/2, weights/2
    phi = np.arange(2*order)*np.pi/order
    mm, pp = np.meshgrid(mu, phi, indexing='ij')
    radial = np.sqrt(1-mm*mm)
    k = np.stack((radial*np.cos(pp), radial*np.sin(pp), mm), axis=-1)
    intensity = np.sum(abs(radiation.amplitude(k))**2, axis=-1)
    result = [float(np.sum(intensity*weight_mu[:, None])*np.pi/order)]
    for na, tilt in zip(apertures, tilts):
        q, w = cap_quadrature(float(na), float(tilt), order)
        result.append(float(w@np.sum(abs(radiation.amplitude(q))**2, axis=-1)) if len(w) else 0.)
    return np.array(result)


def checked_integrals(radiation, apertures, tilts, tolerance=.01):
    previous = integrals(radiation, apertures, tilts, 48)
    for order in (64, 96, 128):
        current = integrals(radiation, apertures, tilts, order)
        changes = abs(current-previous)/np.maximum(abs(current), 1e-300)
        if np.all(changes <= tolerance):
            return current, dict(order=order, relative_changes=changes.tolist(), target_met=True)
        previous = current
    return current, dict(order=order, relative_changes=changes.tolist(), target_met=False)
