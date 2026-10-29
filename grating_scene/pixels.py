"""Exact rectangular-pixel integration of the ideal periodic Jones relay."""
import numpy as np
from scipy.signal import fftconvolve


def integrate_modes(orders, coefficients, config, detector):
    detector.validate()
    x = (np.arange(detector.pixels_x)-(detector.pixels_x-1)/2)*detector.pixel_pitch_um
    y = (np.arange(detector.pixels_y)-(detector.pixels_y-1)/2)*detector.pixel_pitch_um
    if not len(orders):
        return x, y, np.zeros((len(y), len(x)))
    lo, hi = orders.min(axis=0), orders.max(axis=0)
    grid = np.zeros((hi[1]-lo[1]+1, hi[0]-lo[0]+1, 2), complex)
    grid[orders[:, 1]-lo[1], orders[:, 0]-lo[0]] = coefficients
    # Autocorrelation gives coefficients of |Ex|^2+|Ey|^2 without sampling the optical image.
    intensity = sum(fftconvolve(grid[:, :, j], grid[::-1, ::-1, j].conj(), mode='full') for j in (0, 1))
    gx = np.arange(-(hi[0]-lo[0]), hi[0]-lo[0]+1)
    gy = np.arange(-(hi[1]-lo[1]), hi[1]-lo[1]+1)
    lx, ly = config.scene.lx_um*config.magnification, config.scene.ly_um*config.magnification
    active = detector.pixel_pitch_um*np.sqrt(detector.fill_factor)
    intensity *= np.sinc(gx*active/lx)[None, :]*np.sinc(gy*active/ly)[:, None]
    px = np.exp(-2j*np.pi*x[:, None]*gx/lx)
    py = np.exp(-2j*np.pi*y[:, None]*gy/ly)
    integral = py@intensity@px.T*active**2*config.irradiance_w_um2
    scale = max(float(abs(integral).max()), 1e-300)
    if np.max(abs(integral.imag)) > 1e-10*scale or integral.real.min() < -1e-10*scale:
        raise RuntimeError('Periodic pixel integration lost Hermitian positivity')
    return x, y, np.maximum(integral.real, 0.)
