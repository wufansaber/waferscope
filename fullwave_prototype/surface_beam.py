"""Higher-order surface-only ray integration without expanding native FEM expressions."""
import numpy as np

from .beam_inputs import spectrum
from .planar import background_vectors, fresnel


def plane_ray_vectors(p, directions):
    down0, _, _, e0, _, _ = background_vectors(p)
    axis = np.cross(down0, directions)
    incoming = e0+np.cross(axis, e0)+np.cross(axis, np.cross(axis, e0))/(
        1+directions@down0)[:, None]
    rho = np.linalg.norm(directions[:, :2], axis=1)
    s = np.column_stack((-directions[:, 1], directions[:, 0], np.zeros(len(directions))))
    s /= np.where(rho > 1e-14, rho, 1)[:, None]
    s[rho <= 1e-14] = [-1, 0, 0]
    down_p = np.cross(directions, s)
    up = directions*np.array([1, 1, -1])
    rs, rp, _ = fresnel(-directions[:, 2], p.eps_substrate)
    a = np.sum(incoming*s, axis=1)
    b = np.sum(incoming*down_p, axis=1)
    reflected = (a*rs)[:, None]*s+(b*rp)[:, None]*np.cross(up, s)
    return incoming, reflected, np.cross(directions, incoming)


def surface_fields_order(x, y, p, order):
    """Only the numerical integral order changes; physical beam and phase do not."""
    if p.beam_profile != 'gaussian':
        raise ValueError('Surface-only Gaussian integration requires a Gaussian beam')
    directions, weights = spectrum(p, order=order)
    ei, er, hi = plane_ray_vectors(p, directions)
    ei, er, hi = (v*weights[:, None] for v in (ei, er, hi))
    down0 = background_vectors(p)[0]
    central_flux = (1. if p.beam_definition in ('wafer_rectangle', 'measured_wafer') else
                   float(np.real(np.cross(ei.sum(axis=0), hi.sum(axis=0).conj()))@down0))
    if not np.isfinite(central_flux) or central_flux <= 0:
        raise ValueError('Invalid independently integrated central beam flux')
    ei, er, hi = (v/np.sqrt(central_flux) for v in (ei, er, hi))
    groups = np.unique(directions[:, 0], return_index=True)[1]
    groups = np.sort(groups)
    shape = (len(y), len(x), 3)
    total, incident, magnetic = (np.zeros(shape, complex) for _ in range(3))
    for start, stop in zip(groups, np.r_[groups[1:], len(directions)]):
        phase_y = np.exp(1j*p.k0*np.outer(y, directions[start:stop, 1]))
        phase_x = np.exp(1j*p.k0*x*directions[start, 0])[None, :, None]
        total += (phase_y@(ei[start:stop]+er[start:stop]))[:, None, :]*phase_x
        incident += (phase_y@ei[start:stop])[:, None, :]*phase_x
        magnetic += (phase_y@hi[start:stop])[:, None, :]*phase_x
    return total, incident, magnetic


def checked_surface_order(p, width):
    axis = (np.arange(96)-48)*width/96
    previous = surface_fields_order(axis, axis, p, p.beam_order)[0]
    checks, last_order, stable_steps = [], p.beam_order, 0
    for order in (61, 81, 151, 251, 401, 601, 1001, 1601, 2401, 3201):
        if order <= last_order:
            continue
        current = surface_fields_order(axis, axis, p, order)[0]
        change = float(np.linalg.norm(current-previous)/max(np.linalg.norm(current), 1e-300))
        checks.append(dict(orders=[last_order, order], relative_l2_change=change))
        stable_steps = stable_steps+1 if change <= .001 else 0
        if stable_steps >= 2:
            return order, checks, True
        last_order, previous = order, current
    return last_order, checks, False
