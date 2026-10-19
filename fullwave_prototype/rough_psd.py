"""Expected PSD of the actual tapered, piecewise-linear FEM rough patch."""
import numpy as np
from numpy.polynomial.legendre import leggauss

from .rough_surface import _taper_grid, _triangle_indices, taper_effective_area


def height_basis_transform(roughness, frequency, quadrature_order=6):
    """Fourier integrals of nodal hat functions over the projected triangles."""
    frequency = np.asarray(frequency, dtype=float)
    if frequency.ndim != 2 or frequency.shape[1] != 2 or not np.isfinite(frequency).all():
        raise ValueError('Frequency must have shape (directions, 2) and be finite')
    if quadrature_order < 1:
        raise ValueError('Quadrature order must be positive')
    n = roughness.grid_points
    axis = np.linspace(-roughness.width_um/2, roughness.width_um/2, n)
    x, y = np.meshgrid(axis, axis, indexing='ij')
    coordinates = np.column_stack((x.ravel(), y.ravel()))
    indices = _triangle_indices(n)
    corners = coordinates[indices]
    a, b, c = corners[:, 0], corners[:, 1], corners[:, 2]
    determinant = np.abs((b[:, 0]-a[:, 0])*(c[:, 1]-a[:, 1])-
                         (b[:, 1]-a[:, 1])*(c[:, 0]-a[:, 0]))
    nodes, weights = leggauss(quadrature_order)
    nodes, weights = (nodes+1)/2, weights/2
    basis = np.zeros((len(frequency), n*n), dtype=complex)
    for u, wu in zip(nodes, weights):
        for v, wv in zip(nodes, weights):
            points = a+u*(b-a)+(1-u)*v*(c-a)
            weighted_phase = (np.exp(-2j*np.pi*(frequency @ points.T))*
                              (determinant*(1-u)*wu*wv)[None, :])
            for vertex_ids, barycentric in zip(indices.T,
                    ((1-u)*(1-v), u, (1-u)*v)):
                np.add.at(basis.T, vertex_ids, (weighted_phase*barycentric).T)
    return basis


def expected_finite_psd(roughness, theta, phi, wavelength_um, incidence_deg,
                        quadrature_order=6):
    """E[|Fourier(height)|^2]/A_eff, including sampling and edge taper.

    This is a first-order PSD reference, not a second Maxwell solve. It uses
    the same stationary nodal covariance and triangular interpolation as FEM.
    """
    if roughness.model not in ('gaussian', 'abc', 'measured_psd') or roughness.rms_nm <= 0:
        raise ValueError('Finite-patch PSD requires positive statistical roughness')
    if wavelength_um <= 0 or not np.isfinite(wavelength_um):
        raise ValueError('Wavelength must be positive and finite')
    theta, phi = np.asarray(theta), np.asarray(phi)
    polar, azimuth = np.meshgrid(theta, phi, indexing='ij')
    frequency = np.column_stack((np.sin(polar).ravel()*np.cos(azimuth).ravel(),
        np.sin(polar).ravel()*np.sin(azimuth).ravel()-
        np.sin(np.deg2rad(incidence_deg))))/wavelength_um
    basis = height_basis_transform(roughness, frequency, quadrature_order)
    n = roughness.grid_points
    axis = np.linspace(-roughness.width_um/2, roughness.width_um/2, n)
    windowed_basis = basis.reshape((-1, n, n))*_taper_grid(roughness)[None, :, :]
    if roughness.model == 'gaussian':
        kernel = np.exp(-((axis[:, None]-axis[None, :])/
                          roughness.correlation_length_um)**2)
        transformed = np.einsum('ik,fkj->fij', kernel, windowed_basis, optimize=True)
        transformed = np.einsum('fij,jl->fil', transformed, kernel, optimize=True)
        power = np.einsum('fij,fij->f', windowed_basis.conj(), transformed,
                          optimize=True).real
    else:
        from scipy.signal import fftconvolve
        if roughness.model == 'measured_psd':
            from .experimental_surface import psd_covariance_kernel
            kernel = psd_covariance_kernel(roughness)/(roughness.rms_nm/1000)**2
        else:
            from .rough_statistics import abc_correlation
            offsets = np.arange(1-n, n)*roughness.width_um/(n-1)
            dx, dy = np.meshgrid(offsets, offsets, indexing='ij')
            kernel = abc_correlation(np.hypot(dx, dy), roughness.correlation_length_um,
                                     roughness.psd_exponent)
        power = np.empty(len(windowed_basis))
        for start in range(0, len(power), 64):
            block = windowed_basis[start:start+64]
            transformed = fftconvolve(block, kernel[None, :, :],
                                      mode='same', axes=(1, 2))
            power[start:start+64] = np.einsum(
                'fij,fij->f', block.conj(), transformed, optimize=True).real
    power *= (roughness.rms_nm/1000)**2
    return (np.maximum(power, 0)/taper_effective_area(roughness)).reshape(polar.shape)
