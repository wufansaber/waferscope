from dataclasses import replace

import numpy as np
import pytest

from fullwave_prototype.finite_grating import FiniteGrating, trench_vertices
from fullwave_prototype.grating_analysis import relative
from fullwave_prototype.grating_plotting import geometry_faces


@pytest.mark.parametrize('changes', [dict(trench_count=0), dict(trench_count=True),
    dict(trench_count=1.5), dict(top_width_um=.5), dict(bottom_width_um=.6),
    dict(depth_um=.18), dict(length_um=-1), dict(margin_um=0),
    dict(pitch_um=float('nan')), dict(corner_mesh_um=-.01)])
def test_rejects_invalid_grating(changes):
    with pytest.raises(ValueError):
        replace(FiniteGrating(), **changes).validate()


def test_three_trapezoids_have_finite_ends_and_analytic_volume():
    g = FiniteGrating(trench_count=3, bottom_width_um=.15)
    g.validate()
    np.testing.assert_allclose(g.centers, [-.5, 0, .5])
    vertices = trench_vertices(g, 0)
    assert np.ptp(vertices[:, 0]) == pytest.approx(.8)
    assert np.ptp(vertices[:4, 1]) == pytest.approx(.25)
    assert np.ptp(vertices[4:, 1]) == pytest.approx(.15)
    assert g.trench_volume_um3 == pytest.approx(3*.8*.1*.2)
    faces, kinds = geometry_faces(g)
    assert len([k for k in kinds if k == 'groove']) == 15
    assert all(np.isfinite(f).all() for f in faces)


def test_independent_surfaces_are_outside_all_contrast_and_inside_pml():
    g = FiniteGrating()
    for fraction in (.35, .65):
        bounds = np.asarray(g.surface_bounds(fraction))
        assert g.length_um/2 < bounds[0, 1] < g.half_width_um
        assert g.extent_y_um < bounds[1, 1] < g.half_width_um
        assert -g.substrate_depth_um < bounds[2, 0] < -g.depth_um
        assert 0 < bounds[2, 1] < g.air_height_um
    with pytest.raises(ValueError):
        g.surface_bounds(1)


def test_comparison_preserves_absolute_amplitude_and_phase():
    a = np.array([1+1j, 2-3j])
    assert relative(2*a, a) == pytest.approx(1)
    assert relative(1j*a, a) == pytest.approx(np.sqrt(2))
    with pytest.raises(ValueError):
        relative(a, np.ones(3))


def test_saved_run_comparison_rejects_geometry_change(tmp_path):
    import json
    from tools.compare_finite_grating import compare
    base, refined = tmp_path/'base', tmp_path/'refined'
    for folder, geometry in ((base, {'count': 2}), (refined, {'count': 3})):
        folder.mkdir()
        (folder/'report.json').write_text(json.dumps({'geometry': geometry}), encoding='utf-8')
    with pytest.raises(ValueError, match='Identical geometry'):
        compare(base, refined)
