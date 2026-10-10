"""Exact joint Gaussian-process draws for paired height-grid refinement."""
from fractions import Fraction
from functools import lru_cache
import hashlib
import json

import numpy as np


@lru_cache(maxsize=4)
def joint_nodes(levels):
    levels = tuple(sorted(set(levels)))
    if not levels or any(type(n) is not int or not 4 <= n <= 41 for n in levels):
        raise ValueError('Coupled grids require integer levels in [4,41]')
    grids = {n: [(Fraction(i, n-1), Fraction(j, n-1))
                 for i in range(n) for j in range(n)] for n in levels}
    keys = sorted(set(point for grid in grids.values() for point in grid))
    if len(keys) > 2048:
        raise ValueError('Coupled grids exceed 2048 joint nodes; use nested grids such as 9,17,33')
    lookup = {point: index for index, point in enumerate(keys)}
    coordinates = np.array([(float(x)-.5, float(y)-.5) for x, y in keys])
    indices = {n: np.array([lookup[point] for point in grid]) for n, grid in grids.items()}
    return coordinates, indices


@lru_cache(maxsize=2)
def joint_factor(model, levels, width, length, exponent):
    coordinates, _ = joint_nodes(levels)
    from scipy.spatial.distance import cdist
    distance = cdist(coordinates*width, coordinates*width)
    if model == 'gaussian':
        covariance = np.exp(-(distance/length)**2)
    elif model == 'abc':
        from .rough_statistics import abc_correlation
        covariance = abc_correlation(distance, length, exponent)
    else:
        raise ValueError('Paired statistical height refinement supports Gaussian and ABC')
    return np.linalg.cholesky(covariance+1e-12*np.eye(len(coordinates)))


def parent_draw(r):
    levels = tuple(sorted(set(r.sampling_grid_levels)))
    _, indices = joint_nodes(levels)
    if r.grid_points not in indices:
        raise ValueError('Current height grid must belong to the coupled levels')
    factor = joint_factor(r.model, levels, r.width_um, r.correlation_length_um, r.psd_exponent)
    return factor @ np.random.default_rng(r.seed).standard_normal(len(factor))


def height_grid(r):
    values = parent_draw(r)
    _, indices = joint_nodes(tuple(sorted(set(r.sampling_grid_levels))))
    return values[indices[r.grid_points]].reshape((r.grid_points, r.grid_points))


def realization_signature(r):
    if not r.sampling_grid_levels:
        return None
    definition = dict(model=r.model, levels=sorted(set(r.sampling_grid_levels)),
                      width=r.width_um, length=r.correlation_length_um,
                      exponent=r.psd_exponent, seed=r.seed, rms_nm=r.rms_nm)
    digest = hashlib.sha256(json.dumps(definition, sort_keys=True).encode('ascii'))
    digest.update(parent_draw(r).tobytes())
    return digest.hexdigest()
