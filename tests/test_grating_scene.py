from dataclasses import replace

import numpy as np
import pytest

from grating_scene.collection import collect
from grating_scene.config import SceneConfig
from grating_scene.imaging import form_images, mode_coefficients
from grating_scene.near import observation_planes
from rcwa_prototype.comparison import modal_power
from rcwa_prototype.geometry import Scene


def synthetic(scale=1., phase=1.):
    scene = Scene(kind='flat', incidence_deg=0., lx_um=2., harmonics_x=1,
                  harmonics_y=1, grating_method='analytic_li')
    orders = np.array([[0, 0], [1, 0], [-1, 0]])
    kx = orders[:, 0]*scene.wavelength_um/scene.lx_um
    rays = np.column_stack((kx, np.zeros(3), np.sqrt(1-kx*kx)))
    electric = np.zeros((3, 3), complex)
    electric[:, 1] = np.array([.4, .2j, .1+.03j])*scale*phase
    data = dict(orders=orders, directions=rays, electric_reflected_z0=electric,
                propagating=np.ones(3, bool))
    data['reflectance_by_order'] = modal_power(electric, rays, data['propagating'], scene.incidence_deg)
    config = SceneConfig(scene=scene, image_samples=32, numerical_apertures=(.8, .8, 0.),
                         tilts_deg=(0., 0., 0.))
    return config, data


def test_config_roundtrip_and_all_three_near_planes():
    config = SceneConfig()
    assert SceneConfig.from_dict(config.to_dict()) == config
    planes = observation_planes(config)
    assert all(p.shape == (25, 25, 3) for p in planes.values())
    assert np.all(planes['xy'][:, :, 2] == config.near_xy_z_um)
    assert np.all(planes['xz'][:, :, 1] == 0)
    assert np.all(planes['yz'][:, :, 0] == 0)
    assert observation_planes(replace(config, near_enabled=False)) == {}


@pytest.mark.parametrize('changes', [dict(numerical_apertures=(1., .6, .6)),
    dict(tilts_deg=(90., 0., 0.)), dict(image_samples=32), dict(illumination_power_w=-1.)])
def test_invalid_scene_or_unresolved_intensity_grid(changes):
    with pytest.raises(ValueError):
        replace(SceneConfig(), **changes).validate()


def test_three_channels_are_independent_and_union_is_not_double_counted():
    config, ref = synthetic()
    _, sample = synthetic(1.2)
    arrays, report = collect(ref, sample, config)
    a, b, disabled = report['channels']
    assert a['sample_power_w'] == b['sample_power_w']
    assert disabled['sample_power_w'] == 0
    assert report['totals']['sample_union_power_w'] == pytest.approx(a['sample_power_w'])
    assert not arrays['channel_masks'][2].any()


def test_fft_relay_matches_independent_direct_fourier_sum_and_absolute_flux():
    config, ref = synthetic()
    _, sample = synthetic(1.2, phase=np.exp(.2j))
    far, powers = collect(ref, sample, config)
    images, report = form_images(ref, sample, far, config)
    orders, coefficients = mode_coefficients(sample, far['channel_masks'][0], config, 0.)
    xx, yy = np.meshgrid(images['image_x_um'], images['image_y_um'])
    phase = np.exp(-2j*np.pi*(xx[:, :, None]*orders[:, 0]/(config.scene.lx_um*config.magnification)
        + yy[:, :, None]*orders[:, 1]/(config.scene.ly_um*config.magnification)))
    direct = phase@coefficients
    np.testing.assert_allclose(images['ch1_electric_sample'], direct, atol=3e-16)
    assert report['channels'][0]['sample_image_power_w'] == pytest.approx(powers['channels'][0]['sample_power_w'])
    assert report['channels'][0]['parseval_relative_error']['sample'] < 1e-13
    assert np.max(abs(images['ch1_difference_w_um2']-images['ch1_delta_w_um2'])) > 0
    np.testing.assert_allclose(images['ch1_difference_w_um2'],
        images['ch1_delta_w_um2']+images['ch1_interference_w_um2'], atol=1e-20)
    assert not images['ch3_sample_w_um2'].any()


def test_magnification_changes_pixel_density_but_not_cell_power():
    config, data = synthetic()
    far, _ = collect(data, data, config)
    low, a = form_images(data, data, far, config)
    high, b = form_images(data, data, far, replace(config, magnification=40.))
    np.testing.assert_allclose(high['ch1_sample_w_um2']*4, low['ch1_sample_w_um2'], atol=1e-20)
    assert a['channels'][0]['sample_image_power_w'] == pytest.approx(b['channels'][0]['sample_image_power_w'])


def test_tilted_polarization_transport_preserves_order_power():
    config, data = synthetic()
    config = replace(config, tilts_deg=(45., -45., 0.), numerical_apertures=(.95, .95, .95))
    far, powers = collect(data, data, config)
    _, report = form_images(data, data, far, config)
    for p, image in zip(powers['channels'], report['channels']):
        assert image['sample_image_power_w'] == pytest.approx(p['sample_power_w'], rel=1e-13)


def test_zero_particle_contrast_has_zero_image_difference():
    config, data = synthetic()
    far, _ = collect(data, data, config)
    arrays, _ = form_images(data, data, far, config)
    for i in range(1, 4):
        assert np.max(abs(arrays['ch{}_difference_w_um2'.format(i)])) == 0.


def test_relative_phase_only_keeps_intensity_but_not_difference_field():
    config, a = synthetic()
    _, b = synthetic(phase=np.exp(.4j))
    far, _ = collect(a, b, config)
    arrays, _ = form_images(a, b, far, config)
    tolerance = 32*np.finfo(float).eps*arrays['ch1_reference_w_um2'].max()
    np.testing.assert_allclose(arrays['ch1_difference_w_um2'], 0., atol=tolerance)
    assert arrays['ch1_delta_w_um2'].max() > 0.
    np.testing.assert_allclose(arrays['ch1_delta_w_um2'], -arrays['ch1_interference_w_um2'], atol=tolerance)


def test_analyzer_and_transmission_do_not_create_power():
    config, data = synthetic()
    far, _ = collect(data, data, config)
    no_analyzer, _ = form_images(data, data, far, config)
    ax, _ = form_images(data, data, far, replace(config, analyzer='x'))
    ay, _ = form_images(data, data, far, replace(config, analyzer='y'))
    half, _ = form_images(data, data, far, replace(config, transmission=.5))
    np.testing.assert_allclose(ax['ch1_sample_w_um2']+ay['ch1_sample_w_um2'], no_analyzer['ch1_sample_w_um2'], atol=1e-20)
    np.testing.assert_allclose(half['ch1_sample_w_um2']*2, no_analyzer['ch1_sample_w_um2'], atol=1e-20)


@pytest.mark.parametrize('polarization', ['s', 'p'])
def test_symmetric_modal_field_produces_mirrored_independent_channels(polarization):
    config, _ = synthetic()
    config = replace(config, scene=replace(config.scene, harmonics_x=2),
        numerical_apertures=(.8, .8, 0.), tilts_deg=(45., -45., 0.))
    g = np.arange(-2, 3)
    nx = g*config.scene.wavelength_um/config.scene.lx_um
    nz = np.sqrt(1-nx*nx)
    rays = np.column_stack((nx, np.zeros(5), nz))
    polar = np.tile([0., 1., 0.], (5, 1)) if polarization == 's' else np.column_stack((-nz, np.zeros(5), nx))
    coeff = np.array([.12+.05j, .2j, .4, .2j, .12+.05j])
    field = polar*coeff[:, None]
    data = dict(orders=np.column_stack((g, np.zeros(5, int))), directions=rays,
                propagating=np.ones(5, bool), electric_reflected_z0=field)
    data['reflectance_by_order'] = modal_power(field, rays, data['propagating'], 0.)
    far, powers = collect(data, data, config)
    images, _ = form_images(data, data, far, config)
    mirror_index = (-np.arange(config.image_samples)) % config.image_samples
    np.testing.assert_allclose(images['ch1_sample_w_um2'], images['ch2_sample_w_um2'][:, mirror_index], atol=1e-18)
    assert powers['channels'][0]['sample_power_w'] == pytest.approx(powers['channels'][1]['sample_power_w'])
