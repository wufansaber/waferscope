from dataclasses import replace
import json

import numpy as np
import pytest

from grating_scene.config import SceneConfig
from grating_scene.pipeline import DEFAULT_EXECUTABLE, is_within, load_scene, reimage_scene, run_scene
from grating_scene.near import air_modal_closure
from rcwa_prototype.geometry import Scene


@pytest.mark.skipif(not DEFAULT_EXECUTABLE.is_file(), reason='Private S4 stripes runtime required')
def test_complete_scene_retains_native_fields_zero_contrast_and_manifest(tmp_path):
    config = SceneConfig(scene=Scene(harmonics_x=1, harmonics_y=3, slices=2,
        particle_n=1., gap_um=0., grating_method='analytic_li'), near_samples=5, image_samples=32)
    folder = tmp_path/'scene'
    report = run_scene(config, folder)
    loaded, checked, data = load_scene(folder)
    assert loaded == config
    assert checked['validation']['image_parseval_checks_passed']
    assert checked['validation']['air_near_far_closure']['sample']['passed']
    assert not checked['validation']['isolated_defect_certified']
    np.testing.assert_allclose(data['far']['reflected_electric_delta'], 0., atol=3e-12)
    for plane in ('xy', 'xz', 'yz'):
        np.testing.assert_allclose(data['near'][plane+'_delta_electric'], 0., atol=2e-10)
    for name in ('geometry.png', 'near_xy.png', 'near_xz.png', 'near_yz.png',
                 'far_field_2d.png', 'far_field_3d.png', 'images.png'):
        assert (folder/name).stat().st_size > 1000
    assert 'sample/modes.npz' in report['files_sha256']
    with np.load(folder/'sample/modes.npz', allow_pickle=False) as raw:
        phase_test = dict(raw)
    phase_test['near_electric'] *= np.exp(.1j)
    assert not air_modal_closure(phase_test, config)['passed']
    updated = replace(config, numerical_apertures=(0., .4, .8), magnification=40.)
    reused = tmp_path/'reused'
    reimage_scene(updated, folder, reused)
    loaded_updated, checked_updated, derived = load_scene(reused)
    assert loaded_updated == updated
    assert checked_updated['reused_from'] == str(folder.resolve())
    assert not derived['images']['ch1_sample_w_um2'].any()
    assert (reused/'sample/modes.npz').read_bytes() == (folder/'sample/modes.npz').read_bytes()
    with pytest.raises(ValueError, match='require a new S4'):
        reimage_scene(replace(config, scene=replace(config.scene, gap_um=.01)), folder, tmp_path/'bad')
    with pytest.raises(FileExistsError):
        run_scene(config, folder)
    raw = folder/'sample/modes.npz'
    with raw.open('ab') as stream:
        stream.write(b'tampered')
    with pytest.raises(ValueError, match='hash'):
        load_scene(folder)


def test_path_containment_supports_python38_and_rejects_escape(tmp_path):
    assert is_within(tmp_path/'folder/a', tmp_path/'folder')
    assert not is_within(tmp_path/'folder/../secret', tmp_path/'folder')


@pytest.mark.skipif(not DEFAULT_EXECUTABLE.is_file(), reason='Private S4 stripes runtime required')
def test_far_only_control_does_not_invent_near_fields(tmp_path):
    config = SceneConfig(scene=Scene(harmonics_x=1, harmonics_y=3, slices=2,
        grating_method='analytic_li'), near_enabled=False, image_samples=32)
    report = run_scene(config, tmp_path/'far_only')
    _, _, arrays = load_scene(tmp_path/'far_only')
    assert not report['validation']['near_fields_sampled']
    assert arrays['near'] == {}
    assert not (tmp_path/'far_only/near_xy.png').exists()
