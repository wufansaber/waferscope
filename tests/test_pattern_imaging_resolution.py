"""Extended-source imaging must preserve spatial phase and magnification."""
import numpy as np

from fullwave_prototype.imaging import ImagingSettings, image_from_amplitude
from fullwave_prototype.imaging_convergence_audit import dipole_amplitude
from optical_propagation import pupil_geometry
from tools.check_pattern_resolution import translation_check


def test_translated_source_moves_at_actual_magnification():
    result = translation_check()
    np.testing.assert_allclose([r['image_peak_x_um'] for r in result['points']], [-10.64, -21.28])
    assert result['object_space_profile_relative_l2'] < .001


def test_two_separated_coherent_points_are_not_collapsed_to_origin():
    local, k, rotation, mask, du = pupil_geometry(.5, 0, 64)
    base = dipole_amplitude(k[mask], (1, 0, 0))
    amplitude = np.zeros(mask.shape+(2,), complex)
    amplitude[mask] = base*(2*np.cos(2*np.pi/.266*k[mask, 0]*1.064))[:, None]
    settings = ImagingSettings(pupil_samples=64, padding=8, magnification=20)
    result = image_from_amplitude(amplitude, local, k, rotation, du, .266, settings)
    row = result['flux'][len(result['axis_um'])//2]
    peaks = result['axis_um'][np.argsort(row)[-2:]]
    np.testing.assert_allclose(np.sort(peaks), [-21.28, 21.28], atol=1e-10)


def test_large_pattern_view_includes_geometric_arms_without_changing_camera_view():
    from fullwave_prototype.imaging_plotting import image_limits
    from fullwave_prototype.parameters import Parameters
    from dataclasses import asdict
    settings = asdict(ImagingSettings(magnification=20))
    large = asdict(Parameters(defect_shape='cross', defect_length_nm=2000,
                             defect_width_nm=500, gap_nm=10, wavelength_nm=266))
    assert image_limits(settings, large) > 20
    assert image_limits(settings) == 15
    assert image_limits(settings, asdict(Parameters())) == 15
    settings['illumination_fit'] = 'camera'
    assert image_limits(settings, large) == (104, 104)
