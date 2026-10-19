"""Independent ray rotations and direct vector diffraction quadrature (no FFT).

Input is complex (s,p) far-field amplitude in um, with p=k cross s.
Quadrature weights are du*dv, not solid angle; apodization is explicit.
The common axial carrier exp(i*k*z) is omitted, as in the saved image fields.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss


def rotate_between(field, source, target):
    """Rodrigues minimal rotation between non-antiparallel unit directions."""
    cross = np.cross(source, target)
    denominator = 1+np.sum(source*target, axis=-1)
    if np.any(denominator < 1e-12):
        raise ValueError('Antiparallel ray rotation is ambiguous')
    return field+np.cross(cross, field)+np.cross(cross, np.cross(cross, field))/denominator[..., None]


def ray_frame(uv, tilt_deg):
    uv = np.asarray(uv, float)
    if uv.ndim != 2 or uv.shape[1] != 2 or not np.isfinite(uv).all():
        raise ValueError('Finite N by 2 pupil coordinates required')
    radius2 = np.sum(uv**2, axis=1)
    if np.any(radius2 >= 1):
        raise ValueError('Pupil rays must be strictly inside the unit disk')
    t = np.deg2rad(tilt_deg)
    basis = np.array(((np.cos(t), 0, -np.sin(t)), (0, 1, 0), (np.sin(t), 0, np.cos(t))))
    local = np.column_stack((uv, np.sqrt(1-radius2)))
    global_k = local@basis
    phi = np.arctan2(global_k[:, 1], global_k[:, 0])
    s = np.column_stack((-np.sin(phi), np.cos(phi), np.zeros(len(phi))))
    s[np.hypot(global_k[:, 0], global_k[:, 1]) < 1e-12] = (-1, 0, 0)
    return local, global_k, s, np.cross(global_k, s), basis


def polar_cap(na, tilt_deg, radial_order=48, azimuth_count=192):
    """Gauss radial / periodic azimuth quadrature in the projected pupil disk."""
    if not 0 <= na < 1 or radial_order < 1 or azimuth_count < 4:
        raise ValueError('Invalid cap quadrature')
    nodes, weights = leggauss(radial_order)
    radius, dr = (nodes+1)*na/2, weights*na/2
    phi = np.arange(azimuth_count)*2*np.pi/azimuth_count
    uv = np.stack((radius[:, None]*np.cos(phi), radius[:, None]*np.sin(phi)), -1).reshape(-1, 2)
    weight = np.broadcast_to((radius*dr*2*np.pi/azimuth_count)[:, None],
                             (radial_order, azimuth_count)).ravel().copy()
    _, directions, _, _, _ = ray_frame(uv, tilt_deg)
    valid = (directions[:, 2] > 1e-8) & (weight > 0)
    return uv[valid], weight[valid], directions[valid]


def direct_image(amplitude, uv, weights, xy_um, wavelength_um, tilt_deg=0,
                 magnification=20, defocus_um=0, transmission=1, analyzer='none', chunk=128):
    """Sum plane waves with independent Rodrigues rotations and axial E x H*."""
    amplitude, weights = np.asarray(amplitude, complex), np.asarray(weights, float)
    xy = np.asarray(xy_um, float)
    local, direction, s, p, basis = ray_frame(uv, tilt_deg)
    if (amplitude.shape != (len(local), 2) or weights.shape != (len(local),)
            or xy.ndim != 2 or xy.shape[1] != 2):
        raise ValueError('Incompatible amplitude, quadrature or image coordinates')
    if (not all(np.isfinite(a).all() for a in (amplitude, weights, xy)) or np.any(weights < 0)
            or not np.isfinite([wavelength_um, tilt_deg, magnification, defocus_um, transmission]).all()
            or wavelength_um <= 0 or magnification < 1 or not 0 <= transmission <= 1 or chunk < 1):
        raise ValueError('Invalid diffraction input')
    if np.any(direction[:, 2] <= 0):
        raise ValueError('Only the upper hemisphere is supported')
    field = (s*amplitude[:, :1]+p*amplitude[:, 1:])@basis.T
    z = np.broadcast_to((0., 0., 1.), field.shape)
    field = rotate_between(field, local, z)
    field[:, 2] = 0
    if analyzer != 'none':
        angles = {'x': 0, 'y': 90, '45': 45, '135': 135}
        if analyzer not in angles:
            raise ValueError('Unknown analyzer')
        angle = np.deg2rad(angles[analyzer])
        a = np.array((np.cos(angle), np.sin(angle), 0))
        field = (field@a)[:, None]*a
    qxy = -local[:, :2]/magnification
    q = np.column_stack((qxy, np.sqrt(1-np.sum(qxy*qxy, axis=1))))
    e_ray = rotate_between(field, z, q)
    e_ray *= np.sqrt(transmission/(local[:, 2]*q[:, 2]))[:, None]
    e_ray *= np.exp(2j*np.pi/wavelength_um*defocus_um*(q[:, 2]-1))[:, None]
    h_ray = np.cross(q, e_ray)
    weighted = np.column_stack((e_ray, h_ray))*weights[:, None]/(wavelength_um*magnification)
    fields = np.empty((len(xy), 6), complex)
    for start in range(0, len(xy), chunk):
        phase = np.exp(2j*np.pi/wavelength_um*(xy[start:start+chunk]@qxy.T))
        fields[start:start+chunk] = phase@weighted
    e, h = fields[:, :3], fields[:, 3:]
    flux = np.real(np.cross(e, h.conj())[:, 2])
    return e, h, flux
