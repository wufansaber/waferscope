"""Reject unresolved general-envelope illumination before native FEM assembly."""
from dataclasses import replace
import numpy as np
from .background_sampling import electric_background


def check_local_beam(p):
    if p.beam_profile != 'gaussian' or p.beam_definition not in ('wafer_rectangle', 'measured_wafer'):
        return None
    next_order = next((n for n in (31, 41, 61, 81, 101, 151, 251) if n > p.beam_order), None)
    if next_order is None:
        raise ValueError('General-envelope native FEM needs an independent higher-order check')
    x, y, z = np.meshgrid(np.linspace(-.8*p.domain_radius_um, .8*p.domain_radius_um, 3),
                          np.linspace(-.8*p.domain_radius_um, .8*p.domain_radius_um, 3),
                          [-.04, .04, .1], indexing='ij')
    points = np.stack((x, y, z), axis=-1)
    coarse = electric_background(points, p)
    fine = electric_background(points, replace(p, beam_order=next_order))
    change = float(np.linalg.norm(coarse-fine)/max(np.linalg.norm(fine), 1e-300))
    if change > .001:
        raise ValueError('Local coherent beam quadrature changes by {:.3%}; increase beam order '
                         'before native FEM assembly'.format(change))
    return dict(orders=[p.beam_order, next_order], relative_l2_change=change,
                scope='Selected local-domain probes only; not a whole-beam convergence certificate')
