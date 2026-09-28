"""Half-space radiation from volume contrasts on either side of z=0."""
import numpy as np
from numpy.polynomial.legendre import leggauss

from .planar import fresnel


def sphere_quadrature(p, radial_order=10):
    nodes, weights = leggauss(radial_order)
    radius = p.diameter_nm/2000
    r = (nodes+1)*radius/2
    wr = weights*radius/2*r*r
    mu, wm = leggauss(2*radial_order)
    phi = np.arange(4*radial_order)*2*np.pi/(4*radial_order)
    rr, mm, pp = np.meshgrid(r, mu, phi, indexing='ij')
    rho = rr*np.sqrt(1-mm*mm)
    points = np.stack((rho*np.cos(pp), rho*np.sin(pp), rr*mm), axis=-1).reshape(-1, 3)
    points[:, 2] += (p.diameter_nm/2+p.gap_nm)/1000
    w = np.broadcast_to(wr[:, None, None]*wm[None, :, None]*(2*np.pi/len(phi)), rr.shape)
    return points, w.ravel()


def radiation_amplitude(points, polarization_weights, directions, p):
    points = np.asarray(points)
    above = points[:, 2] >= 0
    upper_points, lower_points = points[above], points[~above]
    upper_weights, lower_weights = polarization_weights[above], polarization_weights[~above]
    result = []
    for k in np.asarray(directions):
        rho = np.hypot(k[0], k[1])
        s = np.array((-k[1]/rho, k[0]/rho, 0)) if rho > 1e-12 else np.array((-1., 0, 0))
        kd = np.array((k[0], k[1], -k[2]))
        pup, pdown = np.cross(k, s), np.cross(kd, s)
        rs, rp, q = fresnel(k[2], p.eps_substrate)
        direct = np.exp(-1j*p.k0*(upper_points @ k)) @ upper_weights
        reflected = np.exp(-1j*p.k0*(upper_points @ kd)) @ upper_weights
        # Reflection maps downward s/p amplitudes to upward s/p amplitudes.
        value = np.array((s @ (direct+rs*reflected), pup @ direct+rp*(pdown @ reflected)))
        if len(lower_points):
            # Reciprocal air-incident wave transmitted into the reference substrate.
            ts = (1+rs)*s
            tp = np.asarray((1-rp)*pup, dtype=complex)
            if p.eps_substrate == 1:
                tp = pup.astype(complex)
            else:
                tp[2] = -2*(k[:2] @ pup[:2])/(p.eps_substrate*k[2]+q)
            phase = np.exp(-1j*p.k0*(lower_points @ np.array((k[0], k[1], q))))
            source = phase @ lower_weights
            value += np.array((ts @ source, tp @ source))
        result.append(p.k0**2/(4*np.pi)*value)
    return np.asarray(result)


def angular_grid(step_deg=10):
    theta = np.deg2rad(np.linspace(0, 90, round(90/step_deg)+1))
    phi = np.deg2rad(np.arange(0, 360, step_deg))
    t, f = np.meshgrid(theta, phi, indexing='ij')
    directions = np.stack((np.sin(t)*np.cos(f), np.sin(t)*np.sin(f), np.cos(t)), axis=-1)
    return theta, phi, directions
