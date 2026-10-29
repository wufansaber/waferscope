"""Analytic disk slices for an explicitly selected grcwa Laurent-rule check."""
import math

import numpy as np
from scipy.special import j1

from .geometry import layers


def sphere_stack(scene):
    """Same midpoint slices as S4, without rasterizing the particle boundary."""
    scene.validate()
    if min(scene.lx_um, scene.ly_um) <= scene.diameter_um:
        raise ValueError('Particle must not touch its periodic copies')
    stack, geometry = layers(scene, False)
    radius, dz = scene.diameter_um/2, scene.diameter_um/scene.slices
    volume = 0.
    for i in range(scene.slices):
        r = math.sqrt(radius**2-((i+.5)*dz-radius)**2)
        stack[i] = (dz, dict(radius_um=r) if scene.eps_particle != 1. else None)
        volume += math.pi*r*r*dz
    geometry.update(particle_present=True, voxel_particle_volume_um3=None,
        sliced_particle_volume_um3=volume, exact_particle_volume_um3=4*math.pi*radius**3/3,
        particle_boundary='Analytic disk Fourier coefficients, midpoint vertical slices.')
    return stack, geometry


def disk_convolution(orders, scene, radius_um):
    scene.validate()
    orders = np.asarray(orders)
    if orders.ndim != 2 or orders.shape[1] != 2 or not np.issubdtype(orders.dtype, np.integer):
        raise ValueError('Integer reciprocal orders of shape (N,2) required')
    if not np.isfinite(radius_um) or not 0 < 2*radius_um < min(scene.lx_um, scene.ly_um):
        raise ValueError('Disk must be positive and not touch periodic copies')
    difference = orders[:, None, :]-orders[None, :, :]
    q = 2*np.pi*np.sqrt((difference[..., 0]/scene.lx_um)**2+
                        (difference[..., 1]/scene.ly_um)**2)
    u = q*radius_um
    transform = np.ones_like(u)
    nonzero = u != 0
    transform[nonzero] = 2*j1(u[nonzero])/u[nonzero]
    # grcwa uses an origin at the cell corner; global particle center is (0,0).
    phase = np.exp(-1j*np.pi*np.sum(difference, axis=-1))
    fill = np.pi*radius_um**2/(scene.lx_um*scene.ly_um)
    return np.eye(len(orders), dtype=complex)+(scene.eps_particle-1)*fill*transform*phase


def setup_layer(model, scene, index, radius_um):
    """Only construct epsilon; eigenmodes/scattering remain grcwa operations."""
    from grcwa.rcwa import MakeKPMatrix, SolveLayerEigensystem
    entry = model.id_list[index]
    if entry[0] != 2:
        raise ValueError('Expected an explicit analytic Fourier layer')
    epsilon = disk_convolution(model.G, scene, radius_um)
    epinv = np.linalg.solve(epsilon, np.eye(len(epsilon), dtype=complex))
    zero = np.zeros_like(epsilon)
    ep2 = np.block([[epsilon, zero], [zero, epsilon]])
    kp = MakeKPMatrix(model.omega, 1, epinv, model.kx, model.ky)
    q, phi = SolveLayerEigensystem(model.omega, model.kx, model.ky, kp, ep2)
    if scene.eps_particle.imag == 0 and scene.eps_particle.real > 0:
        q = lossless_branch(q, model.omega)
    model.Patterned_epinv_list[entry[2]] = epinv
    model.Patterned_ep2_list[entry[2]] = ep2
    model.kp_list[index], model.q_list[index], model.phi_list[index] = kp, q, phi


def lossless_branch(q, omega):
    """Do not let roundoff imaginary parts reverse lossless propagating modes."""
    values = np.asarray(q, complex)**2
    scale = max(float(np.max(abs(values))), float(abs(omega)**2))
    near_real = abs(values.imag) <= 256*np.finfo(float).eps*scale
    result = np.asarray(q, complex).copy()
    result[near_real] = np.sqrt(values[near_real].real.astype(complex))
    return result
