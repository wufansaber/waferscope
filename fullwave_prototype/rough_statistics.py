"""Stationary height statistics shared by FEM and its PSD reference."""
from functools import lru_cache
from math import erf, gamma, pi

import numpy as np
from scipy.special import kv


def abc_amplitude(rms_nm, length_um, exponent):
    """ABC amplitude A in um^4 for an infinite 2D PSD with the given RMS."""
    if rms_nm < 0 or length_um <= 0 or exponent <= 2:
        raise ValueError('ABC requires nonnegative RMS, positive B and C > 2')
    return (rms_nm/1000)**2*length_um**2*(exponent-2)/(2*pi)


def abc_correlation(distance_um, length_um, exponent):
    """Unit-variance 2D inverse Fourier transform of NIST's ABC PSD."""
    if length_um <= 0 or exponent <= 2:
        raise ValueError('ABC requires positive B and C > 2')
    distance = np.asarray(distance_um, dtype=float)
    if np.any(distance < 0) or not np.isfinite(distance).all():
        raise ValueError('Distance must be finite and nonnegative')
    order = (exponent-2)/2
    argument = 2*pi*distance/length_um
    result = np.ones_like(argument)
    nonzero = argument > 0
    z = argument[nonzero]
    result[nonzero] = (2**(1-order)/gamma(order))*z**order*kv(order, z)
    return result


@lru_cache(maxsize=4)
def abc_node_covariance(grid_points, width_um, length_um, exponent):
    """Dimensionless covariance of the actual square-grid node heights."""
    axis = np.linspace(-width_um/2, width_um/2, grid_points)
    x, y = np.meshgrid(axis, axis, indexing='ij')
    coordinates = np.column_stack((x.ravel(), y.ravel()))
    distance = np.linalg.norm(coordinates[:, None, :]-coordinates[None, :, :], axis=-1)
    return abc_correlation(distance, length_um, exponent)


@lru_cache(maxsize=4)
def abc_node_factor(grid_points, width_um, length_um, exponent):
    covariance = abc_node_covariance(grid_points, width_um, length_um, exponent)
    return np.linalg.cholesky(covariance+1e-12*np.eye(grid_points**2))


def nyquist_variance_fraction(model, spacing_um, length_um, exponent=3.):
    """Variance within the inscribed Nyquist circle (ABC) or square (Gaussian)."""
    if spacing_um <= 0 or length_um <= 0:
        raise ValueError('Sampling and PSD lengths must be positive')
    nyquist = 1/(2*spacing_um)
    if model == 'gaussian':
        return erf(pi*length_um*nyquist)**2
    if model == 'abc':
        return 1-(1+(length_um*nyquist)**2)**(1-exponent/2)
    raise ValueError('Sampling fraction requires Gaussian or ABC PSD')
