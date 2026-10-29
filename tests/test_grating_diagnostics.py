from dataclasses import replace

import numpy as np
import pytest

from grating_scene.collection import collect
from grating_scene.config import SceneConfig
from grating_scene.diagnostics import scene_diagnostics
from rcwa_prototype.comparison import modal_power
from rcwa_prototype.geometry import Scene


def symmetric(polarization='p'):
    config = SceneConfig(scene=Scene(kind='flat', incidence_deg=0., lx_um=2.,
        harmonics_x=1, harmonics_y=1, grating_method='analytic_li'), polarization=polarization,
        numerical_apertures=(.8, .8, 0.), tilts_deg=(30., -30., 0.), image_samples=32)
    g = np.arange(-1, 2)
    nx = g*config.scene.wavelength_um/config.scene.lx_um
    nz = np.sqrt(1-nx*nx)
    polar = np.tile([0., 1., 0.], (3, 1)) if polarization == 'p' else np.column_stack((nz, np.zeros(3), -nx))
    field = polar*np.array([.2j, .3, .2j])[:, None]
    data = dict(orders=np.column_stack((g, np.zeros(3, int))),
        directions=np.column_stack((nx, np.zeros(3), nz)), propagating=np.ones(3, bool),
        electric_reflected_z0=field)
    data['reflectance_by_order'] = modal_power(field, data['directions'], data['propagating'], 0.)
    far, _ = collect(data, data, config)
    native = {label: dict(reflectance=float(data['reflectance_by_order'].sum()),
        flux_into_substrate=.5, field_vs_library_power_absolute_difference=1e-16)
        for label in ('reference', 'sample')}
    return config, native, far


@pytest.mark.parametrize('polarization', ['s', 'p'])
def test_mirror_diagnostic_obeys_vector_parity_and_does_not_fit(polarization):
    config, native, far = symmetric(polarization)
    audit = scene_diagnostics(far, native, config)
    assert audit['fields']['sample']['mirror_relative_l2'] == 0
    assert audit['fields']['sample']['transverse_relative_l2'] < 1e-15
    assert audit['fields']['sample']['passive_flux_check_passed']
    assert audit['mirror_check_passed']
    assert audit['fields']['delta']['mirror_relative_l2'] is None
    assert audit['side_channel_power_asymmetry']['sample'] == 0
    far['reflected_electric_sample'][0] *= np.exp(.3j)
    changed = scene_diagnostics(far, native, config)
    assert changed['fields']['sample']['mirror_relative_l2'] > .1
    assert not changed['mirror_check_passed']


def test_flux_bookkeeping_rejects_gain_but_does_not_claim_independent_absorption():
    config, native, far = symmetric()
    native['sample']['reflectance'] = .8
    native['sample']['flux_into_substrate'] = .3
    audit = scene_diagnostics(far, native, config)
    assert not audit['fields']['sample']['passive_flux_check_passed']
    assert audit['fields']['sample']['finite_stack_absorption_by_balance'] == pytest.approx(-.1)
    assert 'not independently' in audit['energy_scope']


def test_unmatched_collection_is_not_reported_as_asymmetry_failure():
    config, native, far = symmetric()
    audit = scene_diagnostics(far, native, replace(config, numerical_apertures=(0., .8, 0.)))
    assert not audit['side_channels_mirror_matched']
    assert all(v is None for v in audit['side_channel_power_asymmetry'].values())


def test_longitudinal_contamination_is_detected():
    config, native, far = symmetric()
    far['reflected_electric_sample'] += .1*far['directions']
    audit = scene_diagnostics(far, native, config)
    assert audit['fields']['sample']['transverse_relative_l2'] > .1
