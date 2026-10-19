"""Real Zernike optical-path differences in nm; unit RMS on a full unit disk."""
import numpy as np


def pupil_phase(local, du, wavelength_um, coefficients_nm):
    radius = du*local.shape[0]/2
    x, y = local[..., 0]/radius, local[..., 1]/radius
    r2 = x*x+y*y
    astigmatism, coma, spherical = coefficients_nm
    opd = (astigmatism*np.sqrt(6)*(x*x-y*y)
           + coma*np.sqrt(8)*x*(3*r2-2)
           + spherical*np.sqrt(5)*(6*r2*r2-6*r2+1))
    return np.exp(2j*np.pi*opd/(1000*wavelength_um))
