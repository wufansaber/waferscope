"""Camera footprints, vector-beam fitting and old-source safety."""
from dataclasses import asdict, replace

import numpy as np
import pytest

from detector_model import DetectorParameters
from fullwave_prototype.field_coverage import (camera_footprints, nominal_dimensions,
    illumination_coverage, fit_camera_beam, rough_surface_coverage)
from fullwave_prototype.imaging import ImagingSettings
from fullwave_prototype.imaging_run import build_parser, scattering_settings
from fullwave_prototype.imaging_options import imaging_settings
from fullwave_prototype.parameters import Parameters


DETECTOR = DetectorParameters(pixels_x=16, pixels_y=16)
APERTURES, TILTS = (.5, .5, .5), (45, -45, 0)


def test_camera_footprints_include_projection_and_disabled_channels():
    footprints = camera_footprints(DETECTOR, 20, APERTURES, TILTS)
    assert footprints[0]['width_x_um'] == pytest.approx(5.2*np.sqrt(2))
    assert footprints[0]['height_y_um'] == pytest.approx(5.2)
    assert footprints[1]['width_x_um'] == footprints[0]['width_x_um']
    assert footprints[2]['width_x_um'] == pytest.approx(5.2)
    assert camera_footprints(DETECTOR, 20, (0, 0, .5), TILTS) == footprints[2:]
    assert footprints[0]['wafer_edge_object_defocus_um'] == pytest.approx(2.6)


def test_nominal_gaussian_fit_checks_corners_not_just_axis_edges():
    footprints = camera_footprints(DETECTOR, 20, APERTURES, TILTS)
    length, width = nominal_dimensions(footprints)
    corner = footprints[0]
    ratio = np.exp(-8*((corner['width_x_um']/2/width)**2+(corner['height_y_um']/2/length)**2))
    assert ratio == pytest.approx(.5)


def test_actual_vector_field_is_covered_and_total_power_is_renormalized():
    args = build_parser().parse_args(['--beam-profile', 'gaussian', '--roughness-model', 'gaussian',
        '--rough-grid-points', '9', '--order', '1', '--illumination-fit', 'camera', '--wavelength-nm', '266'])
    p, _ = scattering_settings(args)
    settings = imaging_settings(args, parameters=p)
    assert p.beam_length_um > 12 and p.beam_width_um > 17
    assert settings.irradiance_w_um2*settings.illumination_normalization_area_um2 == pytest.approx(1.)
    assert illumination_coverage(p, DETECTOR, 20, APERTURES, TILTS)['all_channels_covered']
    small = replace(p, beam_length_um=2, beam_width_um=2)
    assert not illumination_coverage(small, DETECTOR, 20, APERTURES, TILTS)['all_channels_covered']
    args.detector_pixels_x = 64
    with pytest.raises(ValueError, match='calculate new scattering'):
        imaging_settings(args, parameters=p)


def test_fit_keeps_beam_algorithm_and_handles_beam_normal_dimensions():
    p = Parameters(beam_profile='gaussian', beam_definition='beam_cross_section',
                   beam_length_um=2, beam_width_um=2)
    fitted, report = fit_camera_beam(p, DETECTOR, 20, APERTURES, TILTS)
    changed = {k for k, v in asdict(p).items() if asdict(fitted)[k] != v}
    assert changed == {'beam_length_um', 'beam_width_um'}
    assert report['all_channels_covered']


def test_surface_taper_must_be_outside_sensor_field():
    settings = ImagingSettings(detector=DETECTOR, illumination_fit='camera')
    assert rough_surface_coverage(settings, dict(width_um=12, taper_width_um=.1))[
        'camera_fields_inside_untapered_surface']
    with pytest.raises(ValueError, match='Rough surface/taper'):
        rough_surface_coverage(settings, dict(width_um=.2, taper_width_um=.01))


def test_legacy_waist_is_not_misreported_as_infinite_plane_illumination():
    report = illumination_coverage(Parameters(beam_waist_um=.1), DETECTOR, 20, APERTURES, TILTS)
    assert not report['all_channels_covered']
    assert report['channels'][0]['minimum_relative_incident_flux'] is None


@pytest.mark.parametrize('target', (0, 1, float('nan')))
def test_invalid_coverage_inputs_fail(target):
    with pytest.raises(ValueError):
        nominal_dimensions(camera_footprints(DETECTOR, 20, APERTURES, TILTS), target)
    with pytest.raises(ValueError):
        nominal_dimensions([], .5)


def test_ui_linked_sizes_and_legacy_load_preserve_physical_beam():
    from PyQt5.QtWidgets import QApplication
    from fullwave_ui.controls import SimulationControls
    app = QApplication.instance() or QApplication([])
    widget = SimulationControls()
    assert widget.fields['illumination-fit'].currentText() == 'camera'
    before = widget.fields['illumination-width-um'].value()
    widget.fields['detector-pixels-x'].setValue(32)
    assert widget.fields['illumination-width-um'].value() == pytest.approx(2*before, abs=.002)
    assert not widget.fields['illumination-width-um'].isEnabled()
    old = asdict(ImagingSettings(detector=DETECTOR, illumination_length_um=2,
                 illumination_width_um=2, illumination_power_w=1))
    widget.load_imaging_settings(old)
    assert widget.fields['illumination-fit'].currentText() == 'manual'
    assert widget.fields['illumination-width-um'].value() == 2
    widget.close()


def test_camera_mode_plots_actual_sensor_roi():
    from fullwave_prototype.imaging_plotting import image_limits
    settings = asdict(ImagingSettings(detector=replace(DETECTOR, pixels_y=8), illumination_fit='camera'))
    assert image_limits(settings) == (52., 26.)
    settings['illumination_fit'] = 'manual'
    assert image_limits(settings) == 15


def test_loaded_physical_inputs_restore_particle_and_stationary_not_local_roughness():
    from PyQt5.QtWidgets import QApplication
    from fullwave_ui.controls import SimulationControls
    app = QApplication.instance() or QApplication([])
    controls = SimulationControls()
    controls.load_physical_parameters(dict(diameter_nm=100, particle_n=1.6, particle_k=.01, order=2),
        dict(model='gaussian', rms_nm=1.2, correlation_length_um=.06, seed=17))
    assert controls.fields['diameter-nm'].value() == 100
    assert controls.fields['particle-k'].value() == .01
    assert controls.fields['rms-nm'].value() == 1.2
    assert controls.fields['roughness-seed'].value() == 17
    assert controls.fields['order'].currentText() == '2'
    controls.close()


def test_wide_field_audit_bounds_only_fft_padding_not_physical_bandwidth():
    from fullwave_prototype.hybrid_surface_audit import audit_padding
    assert audit_padding(1024) == 4
    assert audit_padding(1920) == 2
    assert audit_padding(2048) == 2
    with pytest.raises(ValueError, match='bounded FFT'):
        audit_padding(4096)


def test_wide_gaussian_spectrum_diagnostic_resolves_narrow_angular_support():
    from fullwave_prototype.beam_inputs import truncation_report
    p = Parameters(wavelength_nm=266, incidence_deg=75, beam_profile='gaussian',
        beam_length_um=12.49167305137908, beam_width_um=17.665893445990797)
    report = truncation_report(p)
    assert .999 < report['retained_scalar_spectral_fraction'] <= 1


def test_calculate_routes_camera_coverage_to_whole_surface_hybrid():
    from PyQt5.QtWidgets import QApplication
    from unittest.mock import patch
    from fullwave_ui.widget import FullwaveWidget
    app = QApplication.instance() or QApplication([])
    widget = FullwaveWidget()
    with patch.object(widget, 'start_hybrid_job') as launch:
        widget.start_job()
        launch.assert_called_once()
    with patch.object(widget, 'start_whole_haze') as launch:
        widget.start_ensemble()
        launch.assert_called_once()
    widget.close()
