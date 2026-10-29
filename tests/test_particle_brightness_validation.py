"""Scale-law and geometry controls must not silently change defect type."""
from dataclasses import asdict, replace
import json

import numpy as np
import pytest

from fullwave_prototype.parameters import Parameters
from fullwave_prototype.defects import metadata
from tools.particle_brightness_metrics import local_exponents, geometry_diagnostic
from tools.particle_brightness_phase import constant_phase_relative_l2


def test_rayleigh_power_law_uses_diameter_not_radius_or_linear_intensity():
    d = np.arange(10., 101., 10.)
    np.testing.assert_allclose(local_exponents(d, 3e-20*d**6), 6., atol=1e-12)


def test_invalid_signal_is_not_hidden_by_clipping_or_absolute_values():
    with pytest.raises(ValueError):
        local_exponents([10, 20], [1, -2])
    with pytest.raises(ValueError):
        local_exponents([20, 10], [1, 2])


def test_complex_field_comparison_must_not_fit_intensity_or_position():
    field = np.array([1.+2j, 3.-1j])
    assert constant_phase_relative_l2(field*np.exp(.7j), field) < 1e-14
    assert constant_phase_relative_l2(2*field, field) == pytest.approx(1.)
    assert constant_phase_relative_l2(field*np.exp(1j*np.array([0., 1.])), field) > .1


def test_phase_report_keeps_remaining_field_error_and_power_separate():
    from tools.particle_brightness_note import phase_section
    assert phase_section([]) == []
    text = '\n'.join(phase_section([dict(diameter_nm=100, channel=3,
        old_relative_l2=.2, wafer_reference_relative_l2=.1, pupil_power_relative_change=0.)]))
    assert '| 100 | CH3 | 20 | 10 | 0 |' in text
    assert '相位转换不改变收集功率' in text


def test_background_mesh_refinement_keeps_particle_geometry_and_mesh_fixed():
    from tools.particle_brightness_fem_check import mesh_arguments
    from fullwave_prototype.imaging_run import build_parser, scattering_settings
    values = []
    for scale in (1., .875, .75):
        args = build_parser().parse_args(['--accuracy', 'balanced', '--order', '2',
            '--diameter-nm', '100', '--rms-nm', '0', '--gap-nm', '0']+mesh_arguments(100., scale))
        p, _ = scattering_settings(args)
        assert p.diameter_nm == 100.
        assert p.gap_nm == 0.
        assert p.particle_mesh_um == .01
        values.append(p.substrate_mesh_um)
    np.testing.assert_allclose(values, [.04, .035, .03])
    with pytest.raises(ValueError):
        mesh_arguments(100., float('nan'))


@pytest.mark.parametrize('changed_key,changed_value', [('particle_n', 1.2), ('pml_strength', .5)])
def test_mesh_review_does_not_normalize_power_or_fit_field_amplitude(tmp_path, changed_key, changed_value):
    from tools.particle_brightness_report import mesh_convergence
    for name, scale, amplitude in (('fem_order2_r022_pml2', 1., 1.),
                                   ('fem_order2_r022_pml2_m0875', .875, 2.)):
        folder = tmp_path/name/'d100'
        folder.mkdir(parents=True)
        p = asdict(Parameters(diameter_nm=100, gap_nm=0, air_mesh_um=.08*scale,
                              substrate_mesh_um=.04*scale))
        report = dict(scattering_parameters=p, settings={},
                      channels=[{'with': {'image_power_w': amplitude**2}} for _ in range(3)])
        (folder/'imaging_report.json').write_text(json.dumps(report), encoding='utf-8')
        (folder.parent/'d100.resources.json').write_text(json.dumps(dict(status='completed')), encoding='utf-8')
        np.savez_compressed(folder/'image_fields.npz',
                            **{'ch{}_pupil_with'.format(ch): np.ones((2, 2))*amplitude for ch in (1, 2, 3)})
    rows = mesh_convergence(tmp_path)
    assert len(rows) == 3
    assert all(r['relative_power_change'] == pytest.approx(3.) for r in rows)
    assert all(r['pupil_relative_l2'] == pytest.approx(1.) for r in rows)
    assert not any(r['within_step_targets'] for r in rows)
    report['scattering_parameters'][changed_key] = changed_value
    (folder/'imaging_report.json').write_text(json.dumps(report), encoding='utf-8')
    with pytest.raises(ValueError, match='matching physical'):
        mesh_convergence(tmp_path)


def test_unfinished_fem_run_is_not_counted_as_a_valid_physical_comparison(tmp_path):
    from tools.particle_brightness_report import require_completed
    folder = tmp_path/'d100'
    folder.mkdir()
    for state in (None, 'stopped_memory', 'stopped_timeout', 'failed'):
        if state:
            (tmp_path/'d100.resources.json').write_text(json.dumps(dict(status=state)), encoding='utf-8')
        with pytest.raises(ValueError, match='completed guarded run'):
            require_completed(folder)
    (tmp_path/'d100.resources.json').write_text(json.dumps(dict(status='completed')), encoding='utf-8')
    require_completed(folder)


def test_side_symmetry_compares_absolute_power_and_pixel_peak_separately():
    from tools.particle_brightness_report import side_symmetry
    rows = [dict(incidence_deg=angle, polarization=pol, channel=ch, diameter_nm=d,
        cap_dsc_um2=2., pupil_dsc_um2=2., expected_peak_signal_dn=2. if ch == 1 else 2.02)
        for angle in (0., 75.) for pol in ('p', 's') for ch in (1, 2) for d in (10., 20.)]
    result = side_symmetry(rows)
    assert len(result) == 4
    assert all(r['max_relative_cap_dsc_um2_asymmetry'] == 0 for r in result)
    assert all(r['max_relative_expected_peak_signal_dn_asymmetry'] == pytest.approx(.04/4.02) for r in result)


def test_sphere_diameter_changes_volume_by_cube():
    p = Parameters(defect_shape='sphere')
    small = metadata(replace(p, diameter_nm=10.))
    large = metadata(replace(p, diameter_nm=100.))
    assert large['expected_volume_um3']/small['expected_volume_um3'] == pytest.approx(1000.)


def test_contact_sphere_bottom_remains_zero_when_diameter_changes():
    from fullwave_prototype.defects import center_height_um
    for diameter in range(10, 101, 10):
        p = Parameters(defect_shape='sphere', diameter_nm=float(diameter), gap_nm=0.)
        center = center_height_um(p)
        assert center == pytest.approx(diameter/2000.)
        assert center-diameter/2000. == pytest.approx(0., abs=1e-15)


def test_sphere_diagnostic_overrides_old_gap_without_changing_archived_input():
    p = Parameters(defect_shape='line_x', gap_nm=10.)
    result = geometry_diagnostic(dict(parameters=asdict(p)))
    assert result['sphere_gap_nm'] == 0.
    assert result['sphere_center_heights_nm'] == [5., 25., 50.]
    assert p.gap_nm == 10.


def test_line_dimensions_are_not_changed_by_sphere_diameter_field():
    p = Parameters(defect_shape='line_x', defect_length_nm=200., defect_width_nm=50., defect_height_nm=50.)
    result = geometry_diagnostic(dict(parameters=asdict(p)))
    assert not result['original_geometry_changes_with_diameter']
    assert result['original_geometry'][0]['expected_volume_um3'] == pytest.approx(.0005)
    assert result['sphere_geometry'][0] != result['sphere_geometry'][-1]


@pytest.mark.parametrize('section,key,value', [('parameters', 'substrate_n', 1.),
                                              ('settings', 'transmission', .5)])
def test_comparison_rejects_mismatched_physics_before_loading_arrays(tmp_path, section, key, value):
    from tools.particle_brightness_report import comparisons
    p = asdict(Parameters(wavelength_nm=266, incidence_deg=0, gap_nm=0))
    settings = dict(numerical_apertures=[.5]*3, tilt_degrees=[45., -45., 0.],
                    magnification=20., transmission=1e-6, irradiance_w_um2=.25)
    control = dict(parameters=dict(p), settings=dict(settings))
    (p if section == 'parameters' else settings)[key] = value
    folder = tmp_path/'fem'/'d100'
    folder.mkdir(parents=True)
    (folder/'imaging_report.json').write_text(json.dumps(dict(scattering_parameters=p, settings=settings)),
                                            encoding='utf-8')
    with pytest.raises(ValueError, match='inputs differ'):
        comparisons(tmp_path, dict(rows=[], controls=control))
