"""Exact discrete Gaussian-PSD ensemble for a finite coherent surface field."""
import numpy as np
from scipy.ndimage import map_coordinates

from fullwave_prototype.surface_scattering import radiation_kernel


def covariance_grid(size, spacing, sigma_um, length_um):
    if (not isinstance(size, (int, np.integer)) or isinstance(size, bool) or size < 4
            or not np.isfinite((spacing, sigma_um, length_um)).all()
            or min(spacing, sigma_um, length_um) <= 0):
        raise ValueError('Positive finite covariance grid and PSD parameters required')
    f = np.fft.fftfreq(size, spacing)
    frequency2 = f[:, None]**2+f[None, :]**2
    psd = np.pi*sigma_um**2*length_um**2*np.exp(-(np.pi*length_um)**2*frequency2)
    return np.fft.ifft2(psd).real/spacing**2


def expected_tensor(fields, spacing, sigma_um, length_um):
    """E[A_j(f) A_k(f)*], for the production 2N-periodic parent cropped to N.

    Linear field autocorrelation is multiplied by the height covariance before
    its Fourier transform; this preserves finite-illumination coherence.
    """
    fields = np.asarray(fields)
    if fields.ndim != 3:
        raise ValueError('Finite even square vector field required')
    n = fields.shape[0]
    if (fields.shape != (n, n, 3) or n % 2
            or not np.isfinite((spacing, sigma_um, length_um)).all()
            or min(spacing, sigma_um, length_um) <= 0):
        raise ValueError('Finite even square vector field and positive PSD lengths required')
    if not np.isfinite(fields).all():
        raise ValueError('Non-finite surface field')
    size = 2*n
    covariance = covariance_grid(size, spacing, sigma_um, length_um)
    transforms = [np.fft.fft2(fields[..., j], s=(size, size)) for j in range(3)]
    result = {}
    for j in range(3):
        for k in range(j, 3):
            correlation = np.fft.ifft2(transforms[j]*transforms[k].conj())
            result[j, k] = np.fft.fftshift(np.fft.fft2(covariance*correlation))*spacing**4
    return result


def expected_intensity(tensor, spacing, directions, p):
    q = np.asarray(directions)
    shape = q.shape[:-1]
    flat = q.reshape(-1, 3)
    size = tensor[0, 0].shape[0]
    f = flat[:, :2]/(p.wavelength_nm/1000)
    coordinates = np.array((f[:, 1]*size*spacing+size//2, f[:, 0]*size*spacing+size//2))
    kernel = radiation_kernel(flat, p)
    result = np.zeros(len(flat))
    for (j, k), values in tensor.items():
        moment = map_coordinates(values, coordinates, order=3, prefilter=True,
                                 mode='constant', cval=0)
        factor = np.sum(kernel[..., j]*kernel[..., k].conj(), axis=1)
        result += (factor*moment).real*(1 if j == k else 2)
    if np.min(result) < -1e-8*max(float(result.max()), 1e-300):
        raise ValueError('Covariance reference produced materially negative intensity')
    return np.maximum(result, 0).reshape(shape)
