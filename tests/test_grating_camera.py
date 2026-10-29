from dataclasses import asdict, replace

import numpy as np
import pytest

from detector_model import DetectorParameters, calculate_pixel_detector
from grating_scene.config import SceneConfig
from grating_scene.pixels import integrate_modes
from rcwa_prototype.geometry import Scene


def periodic():
    config = SceneConfig(scene=Scene(kind='flat', incidence_deg=0., lx_um=2., harmonics_x=1,
        harmonics_y=1, grating_method='analytic_li'), image_samples=32)
    orders = np.array([[-1, 0], [0, 0], [1, 1]])
    coeff = np.array([[.1j, .2], [.3, .07j], [.12+.05j, .1j]])
    return config, orders, coeff


def test_analytic_pixel_integral_matches_independent_gauss_quadrature():
    config, orders, coeff = periodic()
    p = DetectorParameters(pixel_pitch_um=3.1, pixels_x=3, pixels_y=2, fill_factor=.63)
    x, y, power = integrate_modes(orders, coeff, config, p)
    nodes, weights = np.polynomial.legendre.leggauss(24)
    active = p.pixel_pitch_um*np.sqrt(p.fill_factor)
    a, b = np.meshgrid(nodes, nodes)
    w = weights[:, None]*weights[None, :]*active**2/4
    for j, yy in enumerate(y):
        for i, xx in enumerate(x):
            phase = np.exp(-2j*np.pi*((xx+a*active/2)[:, :, None]*orders[:, 0]/(config.magnification*config.scene.lx_um)
                +(yy+b*active/2)[:, :, None]*orders[:, 1]/(config.magnification*config.scene.ly_um)))
            field = phase@coeff
            expected = np.sum(w*np.sum(abs(field)**2, axis=-1))*config.irradiance_w_um2
            assert power[j, i] == pytest.approx(expected, rel=2e-14)


def test_one_period_absolute_parseval_and_repeated_camera_roi():
    config, orders, coeff = periodic()
    p = DetectorParameters(pixel_pitch_um=5., pixels_x=8, pixels_y=4)
    _, _, power = integrate_modes(orders, coeff, config, p)
    expected = np.sum(abs(coeff)**2)*config.magnification**2*config.scene.lx_um*config.scene.ly_um*config.irradiance_w_um2
    assert power.sum() == pytest.approx(expected, rel=2e-14)
    _, _, copies = integrate_modes(orders, coeff, config, replace(p, pixels_x=16))
    assert copies.sum() == pytest.approx(2*expected, rel=2e-14)
    np.testing.assert_allclose(copies[:, 4:12], power, rtol=2e-14)
    np.testing.assert_allclose(copies[:, :8], copies[:, 8:], rtol=2e-14)


def test_uniform_pixel_fill_factor_and_disabled_optics():
    config, _, _ = periodic()
    p = DetectorParameters(pixel_pitch_um=5., pixels_x=4, pixels_y=4, fill_factor=.36)
    _, _, power = integrate_modes(np.array([[0, 0]]), np.array([[.2, .3j]]), config, p)
    np.testing.assert_allclose(power, .13*25*.36*config.irradiance_w_um2, rtol=1e-14)
    x, y, zeros = integrate_modes(np.empty((0, 2), int), np.empty((0, 2), complex), config, p)
    assert zeros.shape == (4, 4)
    camera = calculate_pixel_detector(x, y, zeros, 266, p)
    assert camera['signal_e'] == 0
    assert camera['digital_dn'].any()


def test_invalid_integrated_pixel_data_is_rejected():
    p = DetectorParameters(pixels_x=1, pixels_y=1)
    for power in (np.array([[-1.]]), np.array([[np.nan]]), np.zeros((2, 1))):
        with pytest.raises(ValueError):
            calculate_pixel_detector([0], [0], power, 266, p)


def test_native_camera_bundle_and_actual_adc_snr(tmp_path):
    from grating_scene.pipeline import DEFAULT_EXECUTABLE, digest, run_scene
    from grating_scene.camera import load_camera, run_camera
    if not DEFAULT_EXECUTABLE.is_file():
        pytest.skip('Private S4 stripes runtime required')
    from PIL import Image
    config = SceneConfig(scene=Scene(harmonics_x=1, harmonics_y=3, slices=2,
        gap_um=0., grating_method='analytic_li'), near_enabled=False, image_samples=32,
        transmission=1e-6, numerical_apertures=(0., .6, .6))
    source, output = tmp_path/'optics', tmp_path/'camera'
    run_scene(config, source)
    original = digest(source/'sample/modes.npz')
    p = DetectorParameters(enabled=True, pixels_x=4, pixels_y=4, pixel_pitch_um=2.5,
        exposure_ms=.001, gain_e_dn=.5, seed=19, read_noise_e=2.)
    settings = dict(detector=asdict(p), gray_pairs=1000, gray_seed=17)
    report = run_camera(source, output, settings)
    checked, arrays = load_camera(output)
    assert not checked['optical_model_validated']
    assert len(checked['brightness']) == 9
    assert checked['brightness'][0]['adc_max_dn'] == 2**p.bit_depth-1
    assert report['channels'][0]['signal_e'] == 0
    assert digest(output/'optical_scene/sample/modes.npz') == original
    assert digest(source/'sample/modes.npz') == original
    for i in (1, 2, 3):
        prefix = 'ch'+str(i)+'_'
        frames, references = arrays[prefix+'gray_sample_raw_dn'], arrays[prefix+'gray_reference_raw_dn']
        assert frames.dtype == np.uint16 and frames.shape == (1000, 4, 4)
        differences = frames.sum(axis=(1, 2), dtype=float)-references.mean(axis=1).sum(axis=(1, 2))
        channel = checked['channels'][i-1]
        if channel['gray_valid']:
            assert channel['gray_snr'] == pytest.approx(abs(differences.mean())/differences.std(ddof=1))
        with Image.open(output/('ch'+str(i))/'sample_raw.tiff') as image:
            np.testing.assert_array_equal(np.asarray(image), arrays[prefix+'digital_dn'])
    with (output/'camera_data.npz').open('ab') as stream:
        stream.write(b'tampered')
    with pytest.raises(ValueError, match='hash'):
        load_camera(output)
