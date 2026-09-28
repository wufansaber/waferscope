"""Independent homogeneous-sphere Mie oracle, with explicit frame/convention mapping."""
import numpy as np
from scipy.special import spherical_jn, spherical_yn

from .planar import background_vectors


def mie_fields(points, p):
    from miepython.field import eh_near_cartesian
    if p.eps_substrate != 1:
        raise ValueError('Homogeneous Mie oracle requires an air substrate')
    down, _, _, ein, _, _ = background_vectors(p)
    frame = np.column_stack((ein, np.cross(down, ein), down))
    center = np.array((0, 0, (p.diameter_nm/2+p.gap_nm)/1000))
    local = (np.asarray(points)-center) @ frame
    # miepython's lossy-index input is n-ik; its near-field output is exp(-i wt).
    e, h = eh_near_cartesian(p.wavelength_nm/1000, p.diameter_nm/1000,
        complex(p.particle_n, -p.particle_k), 1., *local.T)
    phase = np.exp(1j*p.k0*np.dot(down, center))
    return e.T @ frame.T*phase, h.T @ frame.T*phase


def bessel_mie_cross_sections(p, terms=16):
    """Direct Riccati-Bessel formula, independently of SCATMECH/miepython recurrences."""
    x = p.k0*p.diameter_nm/2000
    m = complex(p.particle_n, p.particle_k)
    n = np.arange(1, terms+1)

    def psi(z):
        return z*spherical_jn(n, z), spherical_jn(n, z)+z*spherical_jn(n, z, True)

    px, dpx = psi(x)
    pm, dpm = psi(m*x)
    xi = px+1j*x*spherical_yn(n, x)
    dxi = dpx+1j*(spherical_yn(n, x)+x*spherical_yn(n, x, True))
    a = (m*pm*dpx-px*dpm)/(m*pm*dxi-xi*dpm)
    b = (pm*dpx-m*px*dpm)/(pm*dxi-m*xi*dpm)
    prefactor = 2*np.pi/p.k0**2
    return dict(csca_um2=float(prefactor*np.sum((2*n+1)*(abs(a)**2+abs(b)**2))),
                cext_um2=float(prefactor*np.sum((2*n+1)*(a+b).real)))
