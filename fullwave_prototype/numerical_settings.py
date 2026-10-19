"""Named starting points, not guarantees of physical accuracy."""
from dataclasses import replace

from .parameters import Parameters


def numerical_parameters(profile='preview', mesh_scale=1., **overrides):
    if profile not in ('preview', 'refined', 'balanced'):
        raise ValueError('Use preview, balanced or refined numerical profile')
    p = Parameters(order=2 if profile == 'preview' else 3,
                   static_condensation=profile != 'preview',
                   pml_strength=1. if profile == 'balanced' else 2.,
                   pml_mesh_scale=1.25 if profile == 'balanced' else 1.,
                   direct_solver='sparsecholesky' if profile == 'balanced' else 'umfpack',
                   air_mesh_um=.08*mesh_scale, substrate_mesh_um=.04*mesh_scale,
                   particle_mesh_um=.015*mesh_scale)
    p = replace(p, **{key: value for key, value in overrides.items() if value is not None})
    p.validate()
    return p
