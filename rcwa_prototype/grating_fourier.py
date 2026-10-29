"""Exact rectangular-line spectrum and Li's 1D inverse-product rule.

This module handles only straight lines; sphere_fourier optionally handles disks.
The grcwa state uses (Hx, Hy) and the rotated electric vector (-Ey, Ex),
so a y-normal interface places the inverse-rule block first in ep2.
"""
from copy import copy

import numpy as np


def line_convolution(orders, scene, reciprocal=False):
    orders = np.asarray(orders)
    if orders.ndim != 2 or orders.shape[1] != 2 or not np.issubdtype(orders.dtype, np.integer):
        raise ValueError('Integer reciprocal orders of shape (N,2) required')
    difference = orders[:, None, :]-orders[None, :, :]
    dy = difference[..., 1]
    allowed = (difference[..., 0] == 0) & (dy % scene.periods_y == 0)
    duty = scene.ridge_width_um/scene.pitch_um
    indicator = np.where(allowed, duty*np.sinc(dy/scene.periods_y*duty)*
                         np.exp(-1j*np.pi*dy), 0.)
    contrast = (1/scene.eps_substrate if reciprocal else scene.eps_substrate)-1
    return np.eye(len(orders), dtype=complex)+contrast*indicator


def material_matrices(orders, scene):
    scene.validate()
    if scene.kind != 'grating':
        raise ValueError('Analytic line matrices require a straight grating')
    epsilon = line_convolution(orders, scene)
    eta = line_convolution(orders, scene, reciprocal=True)
    identity = np.eye(len(orders), dtype=complex)
    epinv = np.linalg.solve(epsilon, identity)
    normal = np.linalg.solve(eta, identity)
    zero = np.zeros_like(epsilon)
    return epinv, np.block([[normal, zero], [zero, epsilon]])


def setup_layer(model, scene, index, block_diagonal=True):
    """Supply material matrices; leave eigenmodes and scattering to grcwa."""
    from grcwa.rcwa import MakeKPMatrix, SolveLayerEigensystem
    entry = model.id_list[index]
    if entry[0] != 2:
        raise ValueError('Expected an explicit analytic Fourier layer')
    orders = np.asarray(model.G)
    n = len(orders)
    epinv = np.zeros((n, n), complex)
    ep2, kp, phi = [np.zeros((2*n, 2*n), complex) for _ in range(3)]
    q = np.zeros(2*n, complex)
    groups = {}
    for i, (gx, gy) in enumerate(orders):
        key = (int(gx), int(gy % scene.periods_y)) if block_diagonal else (0, 0)
        groups.setdefault(key, []).append(i)
    # Straight lines conserve gx and gy modulo the number of primitive periods.
    # Solve those exact independent blocks, not an approximate smaller basis.
    for indices in groups.values():
        take = np.asarray(indices)
        vector = np.r_[take, take+n]
        block_inv, block_ep2 = material_matrices(orders[take], scene)
        block_kp = MakeKPMatrix(model.omega, 1, block_inv, model.kx[take], model.ky[take])
        block_q, block_phi = SolveLayerEigensystem(
            model.omega, model.kx[take], model.ky[take], block_kp, block_ep2)
        epinv[np.ix_(take, take)] = block_inv
        for target, value in ((ep2, block_ep2), (kp, block_kp), (phi, block_phi)):
            target[np.ix_(vector, vector)] = value
        q[vector] = block_q
    model.Patterned_epinv_list[entry[2]] = epinv
    model.Patterned_ep2_list[entry[2]] = ep2
    model.kp_list[index], model.q_list[index], model.phi_list[index] = kp, q, phi
    model._waferscope_exact_line_blocks = block_diagonal


def incident_block(model, scene):
    """Exact reduction for a bare line grating, never for a particle layer."""
    if not getattr(model, '_waferscope_exact_line_blocks', False):
        return model, np.arange(model.nG)
    take = np.flatnonzero((model.G[:, 0] == 0) & (model.G[:, 1] % scene.periods_y == 0))
    vector = np.r_[take, take+model.nG]
    reduced = copy(model)
    reduced.nG = len(take)
    reduced.G, reduced.kx, reduced.ky = model.G[take], model.kx[take], model.ky[take]
    reduced.a0, reduced.bN = model.a0[vector], model.bN[vector]
    for name in ('kp_list', 'phi_list', 'Patterned_ep2_list'):
        setattr(reduced, name, [value[np.ix_(vector, vector)] for value in getattr(model, name)])
    reduced.q_list = [value[vector] for value in model.q_list]
    reduced.Patterned_epinv_list = [value[np.ix_(take, take)] for value in model.Patterned_epinv_list]
    return reduced, take
