"""Exact disk overlap of two finite Fourier fields on matching particle slices."""
import numpy as np

from .contrast_current import disk_transform


def validate_source(source):
    fields = np.asarray(source['source_electric_fourier'])
    wavevectors = np.asarray(source['source_kxy_per_um'])
    z, weight, radius = (np.asarray(source[k]) for k in
        ('source_z_um', 'source_weight_um', 'source_radius_um'))
    if (fields.ndim != 3 or fields.shape[-1] != 3 or wavevectors.shape != (fields.shape[1], 2)
            or z.shape != (fields.shape[0],) or weight.shape != z.shape or radius.shape != z.shape
            or any(not np.all(np.isfinite(a)) for a in (fields, wavevectors, z, weight, radius))
            or np.any(weight <= 0) or np.any(radius <= 0)):
        raise ValueError('Finite Fourier fields, wavevectors and positive particle quadrature required')
    return fields, wavevectors, z, weight, radius


def volume_inner(a, ka, b, kb, weights, radii):
    qx = kb[None, :, 0]-ka[:, None, 0]
    qy = kb[None, :, 1]-ka[:, None, 1]
    total = 0j
    for radius in np.unique(radii):
        take = radii == radius
        disk = disk_transform(qx, qy, radius)
        total += np.einsum('zgk,gh,zhk,z->', a[take].conj(), disk, b[take],
                           weights[take], optimize=True)
    return complex(total)


def compare_source_fields(first, second, epsilon_contrast):
    if not np.isfinite(epsilon_contrast):
        raise ValueError('Finite permittivity contrast required')
    a, ka, za, weights, radii = validate_source(first)
    b, kb, zb, wb, rb = validate_source(second)
    for x, y in ((za, zb), (weights, wb), (radii, rb)):
        if x.shape != y.shape or not np.allclose(x, y, rtol=2e-13, atol=1e-15):
            raise ValueError('Matching physical particle slices and z quadratures required')
    aa = volume_inner(a, ka, a, ka, weights, radii)
    bb = volume_inner(b, kb, b, kb, weights, radii)
    ab = volume_inner(a, ka, b, kb, weights, radii)
    scale = abs(aa)+abs(bb)+2*abs(ab)
    tolerance = 1024*np.finfo(float).eps*scale
    if (aa.real <= 0 or bb.real < 0 or abs(aa.imag) > tolerance or abs(bb.imag) > tolerance):
        raise ValueError('Nonpositive or non-real Fourier energy integral')
    delta = aa.real+bb.real-2*ab.real
    if delta < -tolerance:
        raise ValueError('Negative Fourier difference norm exceeds rounding allowance')

    def moment(fields, k):
        disk = disk_transform(k[None, :, 0], k[None, :, 1], radii[:, None])
        return epsilon_contrast*np.einsum('zgk,zg,z->k', fields, disk, weights, optimize=True)

    ma, mb = moment(a, ka), moment(b, kb)
    norm = float(np.linalg.norm(ma))
    return dict(internal_e_relative_l2=float(np.sqrt(max(0., delta)/aa.real)),
        squared_relative_norm_rounding_allowance=float(tolerance/aa.real),
        first_internal_e_volume_norm_um3=float(aa.real), second_internal_e_volume_norm_um3=float(bb.real),
        first_contrast_moment_um3_real=ma.real.tolist(), first_contrast_moment_um3_imag=ma.imag.tolist(),
        second_contrast_moment_um3_real=mb.real.tolist(), second_contrast_moment_um3_imag=mb.imag.tolist(),
        contrast_moment_relative_change=float(np.linalg.norm(mb-ma)/norm) if norm else None,
        scope='Analytic disk integration of stored Fourier E, existing sliced sphere and z quadrature. '
              'No fitted phase. Relative contrast-current L2 equals internal-E L2 for unchanged nonzero contrast. '
              'Moments are integral((epsilon_particle-1)*E)dV, not collected power or an isolated-limit proof.')
