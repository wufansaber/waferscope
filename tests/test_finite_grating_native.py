from dataclasses import replace

import numpy as np
import pytest

ng = pytest.importorskip('ngsolve')

from fullwave_prototype.finite_grating import FiniteGrating, build_mesh
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.grating_checkpoint import save, load_local
from fullwave_prototype.reciprocity_surface import evaluate_half_field


def small_inputs():
    p = Parameters(domain_radius_um=1.5, order=1, air_mesh_um=.12,
        substrate_mesh_um=.07, pml_thickness_um=.05)
    g = FiniteGrating(trench_count=3, pitch_um=.15, top_width_um=.08,
        bottom_width_um=.05, depth_um=.04, length_um=.3, margin_um=.06,
        air_height_um=.12, substrate_depth_um=.08)
    return p, g


@pytest.mark.parametrize('count', [1, 3])
def test_finite_trapezoid_mesh_volume_and_pmc_boundary(count):
    p, g = small_inputs()
    g = replace(g, trench_count=count)
    mesh = build_mesh(p, g)
    assert 'particle' not in mesh.GetMaterials()
    assert {'outer', 'symmetry'} <= set(mesh.GetBoundaries())
    assert 2*ng.Integrate(1, mesh, definedon=mesh.Materials('pit')) == pytest.approx(g.trench_volume_um3, rel=1e-8)
    physical = sum(float(ng.Integrate(1, mesh, definedon=mesh.Materials(name)))
        for name in ('air', 'substrate', 'pit'))
    assert physical == pytest.approx(2*g.half_width_um**2*(g.air_height_um+g.substrate_depth_um), rel=1e-8)


def test_curved_mesh_checkpoint_preserves_electric_and_curl_in_physical_domain(tmp_path):
    p, g = small_inputs()
    mesh = build_mesh(p, g)
    space = ng.HCurl(mesh, order=1, complex=True)
    field = ng.GridFunction(space)
    field.Set(ng.CF((ng.y+1j*ng.z, ng.x-2j*ng.z, 1+ng.x*ng.y)))
    points = np.array([[.02, .03, .06], [-.04, .02, -.06], [.06, .15, -.02]])
    expected = evaluate_half_field(mesh, field, points, p)
    save(tmp_path, mesh, field, p, g, {'algebraic_residual': 0})
    restored_mesh, restored_field, _ = load_local(tmp_path)
    actual = evaluate_half_field(restored_mesh, restored_field, points, p)
    for a, b in zip(actual, expected):
        np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-12)
    path = tmp_path/'field.ngs'
    path.write_bytes(path.read_bytes()+b'x')
    with pytest.raises(ValueError, match='hash mismatch'):
        load_local(tmp_path)
