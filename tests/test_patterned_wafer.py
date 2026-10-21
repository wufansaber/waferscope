from dataclasses import replace

import numpy as np
import pytest

from fullwave_prototype.parameters import Parameters
from fullwave_prototype.patterned_wafer import LinePatch, surface_triangles, mirror_sources


def test_finite_line_dimensions_and_real_vertical_walls():
    patch = LinePatch()
    patch.validate(Parameters())
    assert patch.trench_volume_um3 == pytest.approx(.14)
    triangles = surface_triangles(patch)
    assert triangles.min(axis=(0, 1))[2] == pytest.approx(-.1)
    assert triangles.max(axis=(0, 1))[2] == 0
    areas = np.linalg.norm(np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0]), axis=1)/2
    assert np.all(areas > 0)
    assert any(np.ptp(t[:, 2]) > 0 for t in triangles)


@pytest.mark.parametrize('change', [dict(pitch_um=2), dict(trench_depth_um=.2),
    dict(trench_width_um=1), dict(line_length_um=3), dict(trench_depth_um=float('nan'))])
def test_pattern_rejects_clipping_or_invalid_dimensions(change):
    with pytest.raises(ValueError):
        replace(LinePatch(), **change).validate(Parameters())


def test_pmc_reconstruction_reflects_vector_not_only_positions():
    points = np.array([[.1, .2, -.03]])
    fields = np.array([[1+2j, 3+4j, 5+6j]])
    full_points, full_fields = mirror_sources(points, fields)
    np.testing.assert_array_equal(full_points[1], [-.1, .2, -.03])
    np.testing.assert_array_equal(full_fields[1], [-1-2j, 3+4j, 5+6j])


@pytest.mark.parametrize('diameter', [100, 200])
def test_requested_40x_500nm_case_applies_to_physics_and_optics(diameter):
    from tools.patterned_wafer_case import build_parser, case_inputs
    args = build_parser().parse_args(['--diameter-nm', str(diameter), '--pitch-nm', '500',
        '--line-width-nm', '250', '--magnification', '40', '--na', '.6',
        '--line-length-um', '.8', '--half-width-um', '.55', '--air-height-um', '.35',
        '--exposure-ms', '.001', '--source-order', '7'])
    p, patch, optics = case_inputs(args)
    assert p.diameter_nm == diameter
    assert patch.pitch_um == .5 and patch.trench_width_um == .25
    assert patch.pitch_um-patch.trench_width_um == .25
    assert optics.magnification == 40 and optics.numerical_apertures == (.6,)*3
    assert optics.focal_length_mm*optics.magnification == 800
    assert p.beam_profile == 'plane' and patch.trench_depth_um == .1
    assert optics.detector.exposure_ms == .001


def test_requested_particle_must_fit_ridge_before_starting_solver():
    from tools.patterned_wafer_case import build_parser, case_inputs
    args = build_parser().parse_args(['--diameter-nm', '300', '--pitch-nm', '500',
                                     '--line-width-nm', '250', '--air-height-um', '.5'])
    with pytest.raises(ValueError, match='central ridge'):
        case_inputs(args)


def test_patterned_case_uses_explicit_high_order_integration():
    from tools.patterned_wafer_case import build_parser, case_inputs
    parser = build_parser()
    p, _, _ = case_inputs(parser.parse_args([]))
    assert p.integration_bonus == 10
    p, _, _ = case_inputs(parser.parse_args(['--order', '2', '--pml-order', '1',
                                            '--integration-bonus', '10']))
    assert p.order == 2 and p.pml_order == 1 and p.static_condensation
    assert p.integration_bonus == 10
