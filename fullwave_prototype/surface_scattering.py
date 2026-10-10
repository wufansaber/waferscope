"""First-order rough-interface radiation, with coherent vector illumination."""
import numpy as np
from scipy.ndimage import map_coordinates, spline_filter

from .planar import fresnel


def radiation_kernel(directions, p):
    """Surface-displacement derivative: parallel E and continuous normal D."""
    q = np.asarray(directions, dtype=float)
    phi = np.arctan2(q[..., 1], q[..., 0])
    s = np.stack((-np.sin(phi), np.cos(phi), np.zeros_like(phi)), axis=-1)
    s = np.where((np.hypot(q[..., 0], q[..., 1]) < 1e-12)[..., None], (-1., 0., 0.), s)
    down = q.copy()
    down[..., 2] *= -1
    rs, rp, _ = fresnel(q[..., 2], p.eps_substrate)
    ss = (1+rs)[..., None]*s
    pp = np.cross(q, s)+rp[..., None]*np.cross(down, s)
    # A thin displaced interface has E_normal = D_normal / eps_substrate.
    pp[..., 2] /= p.eps_substrate
    return np.stack((ss, pp), axis=-2)*p.k0**2*(p.eps_substrate-1)/(4*np.pi)


def surface_fields(x, y, p):
    """Exact ray quadrature on a tensor grid, factoring common X phases."""
    from .illumination import illumination_rays
    groups = {}
    for down, _, _, incoming, reflected, _ in illumination_rays(p):
        groups.setdefault(float(down[0]), []).append((down, incoming, reflected))
    shape = (len(y), len(x), 3)
    total, incident, magnetic = (np.zeros(shape, complex) for _ in range(3))
    for qx, group in groups.items():
        down = np.array([item[0] for item in group])
        ei = np.array([item[1] for item in group])
        er = np.array([item[2] for item in group])
        phase_y = np.exp(1j*p.k0*np.outer(y, down[:, 1]))
        phase_x = np.exp(1j*p.k0*x*qx)[None, :, None]
        total += (phase_y@(ei+er))[:, None, :]*phase_x
        incident += (phase_y@ei)[:, None, :]*phase_x
        magnetic += (phase_y@np.cross(down, ei))[:, None, :]*phase_x
    return total, incident, magnetic


class SurfaceRadiation:
    def __init__(self, x, y, weighted_fields, p, padding=4):
        self.p = p
        self.dx, self.dy = x[1]-x[0], y[1]-y[0]
        ny, nx = weighted_fields.shape[:2]
        if nx % 2 or ny % 2 or padding not in (2, 4):
            raise ValueError('Use even tensor grids and surface FFT padding 2 or 4')
        if (weighted_fields.shape != (len(y), len(x), 3) or
                not np.isfinite(weighted_fields).all() or self.dx <= 0 or self.dy <= 0 or
                not np.allclose(np.diff(x), self.dx) or not np.allclose(np.diff(y), self.dy) or
                not np.isclose(x[nx//2], 0) or not np.isclose(y[ny//2], 0)):
            raise ValueError('Surface Fourier fields require finite, uniform, origin-centered XY grids')
        self.nx, self.ny = nx*padding, ny*padding
        self.spectra = []
        for component in range(3):
            values = np.pad(weighted_fields[..., component],
                ((self.ny//2-ny//2, self.ny//2-ny//2),
                 (self.nx//2-nx//2, self.nx//2-nx//2)))
            transformed = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(values)))
            self.spectra.append(spline_filter(transformed*self.dx*self.dy, order=3,
                                             output=np.complex128))

    def amplitude(self, directions):
        directions = np.asarray(directions)
        shape = directions.shape[:-1]
        frequency = directions.reshape((-1, 3))[:, :2]/(self.p.wavelength_nm/1000)
        coordinates = np.array((frequency[:, 1]*self.ny*self.dy+self.ny//2,
                                frequency[:, 0]*self.nx*self.dx+self.nx//2))
        transformed = np.column_stack([map_coordinates(a, coordinates, order=3,
            mode='constant', cval=0, prefilter=False) for a in self.spectra])
        kernel = radiation_kernel(directions.reshape((-1, 3)), self.p)
        return np.einsum('nij,nj->ni', kernel, transformed).reshape(shape+(2,))


def direct_amplitude(x, y, weighted_fields, directions, p):
    """Independent spatial Fourier sum for selected directions, without FFT interpolation."""
    transforms = []
    for q in directions:
        phase = np.exp(-1j*p.k0*y*q[1])[:, None]*np.exp(-1j*p.k0*x*q[0])[None, :]
        transforms.append(np.sum(weighted_fields*phase[..., None], axis=(0, 1))*(x[1]-x[0])*(y[1]-y[0]))
    return np.einsum('nij,nj->ni', radiation_kernel(directions, p), transforms)
