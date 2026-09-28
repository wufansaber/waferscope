"""Air/half-space Fresnel background. exp(-i omega t), incoming from -KY."""
import numpy as np


def fresnel(cosine, eps):
    cosine = np.asarray(cosine)
    # Avoid losing cos(theta)^2 to cancellation at grazing incidence in air.
    q = np.sqrt((eps-1) + cosine**2 + 0j)
    q = np.where((q.imag < 0) | ((q.imag == 0) & (q.real < 0)), -q, q)
    if eps == 1:
        return np.zeros_like(q), np.zeros_like(q), q
    rs = (cosine-q)/(cosine+q)
    rp = (eps*cosine-q)/(eps*cosine+q)
    return rs, rp, q


def background_vectors(parameters):
    t = np.deg2rad(parameters.incidence_deg)
    s, c = np.sin(t), np.cos(t)
    rs, rp, q = fresnel(c, parameters.eps_substrate)
    down, up = np.array((0, s, -c)), np.array((0, s, c))
    kt = np.array((0, s, -q), dtype=complex)
    if parameters.polarization == 's':
        ein = np.array((-1., 0, 0))
        er = rs*ein
        et = (1+rs)*ein
    else:
        ein = np.array((0., c, s))
        er = rp*np.array((0., -c, s))
        # Tangential E continuity and k dot E = 0 in the substrate.
        et = np.array((0, c*(1-rp), 2*s*c/(parameters.eps_substrate*c+q)), dtype=complex)
    return down, up, kt, ein, er, et


def background_fields(points, parameters):
    down, up, kt, ein, er, et = background_vectors(parameters)
    points = np.asarray(points)
    pi = np.exp(1j*parameters.k0*(points @ down))
    pr = np.exp(1j*parameters.k0*(points @ up))
    pt = np.exp(1j*parameters.k0*(points @ kt))
    eair = pi[..., None]*ein + pr[..., None]*er
    hair = pi[..., None]*np.cross(down, ein) + pr[..., None]*np.cross(up, er)
    esub = pt[..., None]*et
    hsub = pt[..., None]*np.cross(kt, et)
    return (np.where(points[..., 2, None] >= 0, eair, esub),
            np.where(points[..., 2, None] >= 0, hair, hsub))
