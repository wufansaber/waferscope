"""Defect topology, dimensions, archive/UI round trips and native CAD volumes."""
from dataclasses import asdict, replace
import os

import numpy as np
import pytest

from fullwave_prototype.defects import (boundary_distance_um, cad_solid, center_height_um,
    footprint_nm, metadata, required_radius_um, section_nm)
from fullwave_prototype.parameters import Parameters


@pytest.mark.parametrize('shape,volume', [('line_x', .0005), ('line_y', .0005), ('cross', .000875)])
def test_exact_union_volume_and_bounds(shape, volume):
    p = Parameters(defect_shape=shape, gap_nm=10)
    p.validate()
    assert metadata(p)['expected_volume_um3'] == pytest.approx(volume)
    polygon = footprint_nm(p)
    area = .5*abs(np.sum(np.cross(polygon, np.roll(polygon, -1, axis=0))))
    assert area*p.defect_height_nm/1e9 == pytest.approx(volume)
    assert center_height_um(p) == pytest.approx(.035)
    assert required_radius_um(p) < p.domain_radius_um


def test_rotation_and_square_cross():
    horizontal = Parameters(defect_shape='line_x', defect_rotation_deg=90)
    vertical = replace(horizontal, defect_shape='line_y', defect_rotation_deg=0)
    np.testing.assert_allclose(footprint_nm(horizontal), footprint_nm(vertical), atol=1e-12)
    cross = replace(horizontal, defect_shape='cross', defect_length_nm=50)
    assert len(footprint_nm(cross)) == 4
    assert metadata(cross)['expected_volume_um3'] == pytest.approx(.000125)


@pytest.mark.parametrize('changes', [dict(defect_shape='triangle'), dict(defect_height_nm=0),
    dict(defect_shape='cross', defect_width_nm=201), dict(defect_rotation_deg=float('nan')),
    dict(defect_shape='line_x', defect_length_nm=1000)])
def test_reject_invalid_dimensions_and_clipped_domain(changes):
    with pytest.raises(ValueError):
        replace(Parameters(), **changes).validate()


def test_distances_and_cut_planes_use_real_shape():
    p = Parameters(defect_shape='cross', gap_nm=10)
    points = np.array([[0, 0, .035], [.08, .08, .035], [.1, 0, .035], [0, 0, .07]])
    np.testing.assert_allclose(boundary_distance_um(points, p), [.025, .055, 0, .01], atol=1e-14)
    np.testing.assert_allclose(section_nm(p, 'xz'), [[-100, 10], [100, 10], [100, 60], [-100, 60]])
    assert len(section_nm(p, 'xy', 35)) == 12
    assert section_nm(p, 'xy', 70) is None
    from fullwave_prototype.near_planes import make_planes
    xy = next(q for q in make_planes(p) if q.name == 'xy')
    assert xy.offset == pytest.approx(.035)
    for plane in make_planes(p):
        assert np.max(np.linalg.norm(plane.points, axis=1)) < p.domain_radius_um


def test_spherical_references_do_not_accept_patterns():
    from fullwave_prototype.reference import calculate
    from fullwave_prototype.far_field import sphere_quadrature
    p = Parameters(defect_shape='cross')
    for call in (lambda: calculate(p, [0], [0]), lambda: sphere_quadrature(p)):
        with pytest.raises(ValueError, match='sphere'):
            call()


def test_long_line_xy_window_and_corner_bounds():
    from fullwave_prototype.near_planes import make_planes
    p = Parameters(defect_shape='line_x', defect_length_nm=1000, domain_radius_um=.6, gap_nm=10)
    xy = next(q for q in make_planes(p) if q.name == 'xy')
    assert max(xy.u) > .5
    assert max(xy.v) > .025
    assert np.max(np.linalg.norm(xy.points, axis=1)) < p.domain_radius_um


def test_ui_cli_and_saved_geometry_roundtrip():
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt5.QtWidgets import QApplication
    from fullwave_ui.controls import SimulationControls
    from fullwave_prototype.imaging_run import build_parser, scattering_settings
    app = QApplication.instance() or QApplication([])
    controls = SimulationControls()
    f = controls.fields
    f['defect-shape'].setCurrentText('cross')
    f['defect-length-nm'].setValue(240)
    f['defect-width-nm'].setValue(40)
    f['defect-height-nm'].setValue(35)
    f['defect-rotation-deg'].setValue(25)
    assert not f['diameter-nm'].isEnabled()
    assert f['defect-height-nm'].isEnabled()
    args = build_parser().parse_args(controls.arguments())
    p, _ = scattering_settings(args)
    assert (p.defect_shape, p.defect_length_nm, p.defect_width_nm,
            p.defect_height_nm, p.defect_rotation_deg) == ('cross', 240, 40, 35, 25)
    controls.load_physical_parameters(asdict(p))
    assert f['defect-rotation-deg'].value() == 25
    controls.load_physical_parameters({'diameter_nm': 100})
    assert f['defect-shape'].currentText() == 'sphere'
    assert f['diameter-nm'].isEnabled()
    assert not f['defect-height-nm'].isEnabled()
    assert '--defect-shape' not in controls.arguments(imaging_only=True)
    controls.close()
    app.processEvents()


@pytest.mark.parametrize('name,shape', [('line X', 'line_x'), ('line Y', 'line_y'), ('cross', 'cross')])
def test_pattern_presets_have_bounded_preview_settings(name, shape):
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt5.QtWidgets import QApplication
    from fullwave_ui.controls import SimulationControls
    from fullwave_prototype.imaging_run import build_parser, scattering_settings
    app = QApplication.instance() or QApplication([])
    controls = SimulationControls()
    controls.preset.setCurrentText('266 nm '+name+' preview')
    args = build_parser().parse_args(controls.arguments())
    p, roughness = scattering_settings(args)
    assert p.defect_shape == shape
    assert p.order == 1 and p.particle_mesh_um == .015
    assert p.beam_profile == 'plane' and roughness.rms_nm == 0
    assert args.illumination_fit == 'manual' and args.transmission == .000001
    controls.close()
    app.processEvents()


@pytest.mark.parametrize('shape', ['line_x', 'line_y', 'cross'])
@pytest.mark.parametrize('gap', [0., 10.])
def test_native_cad_volume(shape, gap):
    pytest.importorskip('netgen.occ')
    p = Parameters(defect_shape=shape, defect_rotation_deg=17, gap_nm=gap)
    solid = cad_solid(p, 1000)
    assert solid.mass/1e9 == pytest.approx(metadata(p)['expected_volume_um3'], rel=1e-12)


@pytest.mark.parametrize('shape', ['line_x', 'line_y', 'cross'])
@pytest.mark.parametrize('gap', [0., 10.])
def test_native_contact_pattern_mesh_material_volume(shape, gap):
    ng = pytest.importorskip('ngsolve')
    from fullwave_prototype.solver import build_mesh
    p = Parameters(defect_shape=shape, gap_nm=gap, order=1, particle_mesh_um=.02,
                   air_mesh_um=.12, substrate_mesh_um=.06)
    mesh = build_mesh(p)
    volume = ng.Integrate(1, mesh, definedon=mesh.Materials('particle'))
    assert volume == pytest.approx(metadata(p)['expected_volume_um3'], rel=1e-9)
