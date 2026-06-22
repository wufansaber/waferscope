"""Zernike polynomial evaluation using Noll-style Z1...Z36 ordering."""

from math import factorial

import numpy as np


def noll_index_to_nm(index):
    """Map a one-based Noll index to radial order n and azimuthal order m."""
    if index < 1:
        raise ValueError("Zernike index must be one-based and positive")

    remaining = index
    radial_order = 0
    while remaining > radial_order + 1:
        remaining -= radial_order + 1
        radial_order += 1

    if radial_order % 2 == 0:
        orders = [0]
        for magnitude in range(2, radial_order + 1, 2):
            orders.extend((-magnitude, magnitude))
    else:
        orders = []
        for magnitude in range(1, radial_order + 1, 2):
            orders.extend((-magnitude, magnitude))

    return radial_order, orders[remaining - 1]


def radial_polynomial(radial_order, azimuthal_order, radius):
    """Evaluate the radial Zernike polynomial R_n^|m|."""
    magnitude = abs(azimuthal_order)
    if (radial_order - magnitude) % 2:
        return np.zeros_like(radius)

    result = np.zeros_like(radius, dtype=float)
    upper = (radial_order - magnitude) // 2
    for summation_index in range(upper + 1):
        coefficient = (
            (-1) ** summation_index
            * factorial(radial_order - summation_index)
            / (
                factorial(summation_index)
                * factorial((radial_order + magnitude) // 2 - summation_index)
                * factorial((radial_order - magnitude) // 2 - summation_index)
            )
        )
        result += coefficient * radius ** (radial_order - 2 * summation_index)
    return result


def zernike_mode(index, radius, angle):
    """Evaluate a normalized real-valued Zernike mode on a polar grid."""
    radial_order, azimuthal_order = noll_index_to_nm(index)
    radial = radial_polynomial(radial_order, azimuthal_order, radius)

    if azimuthal_order == 0:
        return np.sqrt(radial_order + 1.0) * radial
    normalization = np.sqrt(2.0 * (radial_order + 1.0))
    if azimuthal_order > 0:
        return normalization * radial * np.cos(azimuthal_order * angle)
    return normalization * radial * np.sin(abs(azimuthal_order) * angle)


def evaluate_zernike_wavefront(coefficients, grid_size=181):
    """Return x, y and wavefront arrays for coefficients Z1...ZN in waves."""
    coordinate = np.linspace(-1.0, 1.0, grid_size)
    x, y = np.meshgrid(coordinate, coordinate)
    radius = np.hypot(x, y)
    angle = np.arctan2(y, x)
    pupil = radius <= 1.0
    wavefront = np.zeros_like(radius)

    for index, coefficient in enumerate(coefficients, start=1):
        if coefficient:
            wavefront += float(coefficient) * zernike_mode(index, radius, angle)

    wavefront[~pupil] = np.nan
    valid = wavefront[pupil]
    metrics = {
        "pv": float(np.max(valid) - np.min(valid)),
        "rms": float(np.std(valid)),
        "minimum": float(np.min(valid)),
        "maximum": float(np.max(valid)),
    }
    return x, y, wavefront, metrics
