"""CAD/physical-domain checks; these do not replace Maxwell convergence."""
from dataclasses import replace

import numpy as np
import pytest

ng = pytest.importorskip('ngsolve')

from fullwave_prototype.parameters import Parameters
from fullwave_prototype.patterned_wafer import LinePatch, build_mesh


@pytest.mark.parametrize('patterned', [True, False])
def test_shared_half_domain_preserves_volume_and_excludes_pml_from_source(patterned):
    p = Parameters(wavelength_nm=266, diameter_nm=100, gap_nm=10, domain_radius_um=1.5,
        order=1, air_mesh_um=.15, substrate_mesh_um=.07, particle_mesh_um=.02,
        pml_thickness_um=.08)
    patch = LinePatch(pitch_um=.5, trench_width_um=.25, line_length_um=.8,
                      half_width_um=.55, air_height_um=.35)
    mesh = build_mesh(p, patch, half_domain=True, patterned=patterned)
    names = set(mesh.GetMaterials())
    assert ('pit' in names) == patterned
    assert {'air', 'substrate', 'particle', 'pml_air', 'pml_substrate'} <= names
    volumes = {name: float(ng.Integrate(1, mesh, definedon=mesh.Materials(name))) for name in names}
    physical = sum(v for name, v in volumes.items() if not name.startswith('pml_'))
    expected = 2*patch.half_width_um**2*(patch.air_height_um+patch.substrate_depth_um)
    assert physical == pytest.approx(expected, rel=2e-6)
    assert volumes['particle']*2 == pytest.approx(4*np.pi*(.05**3)/3, rel=1e-4)
    if patterned:
        assert volumes['pit']*2 == pytest.approx(patch.trench_volume_um3, rel=1e-9)
    assert 'symmetry' in mesh.GetBoundaries()
    assert 'outer' in mesh.GetBoundaries()
