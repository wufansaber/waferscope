import json

import pytest

from grating_scene.pipeline import digest, write_json
from grating_scene.resources import MODEL, calibrated_estimate, native_estimate


def profile(tmp_path):
    executable = tmp_path/'native.exe'
    executable.write_bytes(b'test fixture only; never executed')
    rows = [dict(modes=n, layers=l, image_samples=256,
                 peak_accounted_gib=native_estimate(n, l)['planning_gib']/1.4)
            for n, l in ((231, 5), (441, 17), (873, 11))]
    path = tmp_path/'profile.json'
    write_json(path, dict(model=MODEL, executable_sha256=digest(executable),
        status='calibration_checks_passed', observations=rows))
    return executable, path


def test_scoped_native_plan_requires_matching_build_and_retains_profile_hash(tmp_path):
    exe, path = profile(tmp_path)
    result = calibrated_estimate(path, exe, 715, 17, 256)
    assert result['planning_gib'] == pytest.approx(3.425252914428711)
    assert result['profile_sha256'] == digest(path)
    exe.write_bytes(b'changed binary')
    with pytest.raises(ValueError, match='executable'):
        calibrated_estimate(path, exe, 715, 17, 256)


def test_calibration_without_headroom_is_rejected(tmp_path):
    exe, path = profile(tmp_path)
    report = json.loads(path.read_text())
    row = report['observations'][0]
    row['peak_accounted_gib'] = native_estimate(row['modes'], row['layers'])['planning_gib']
    write_json(path, report)
    with pytest.raises(ValueError, match='headroom'):
        calibrated_estimate(path, exe, 715, 17, 256)


def test_calibration_cannot_be_extrapolated_to_arbitrary_dimensions(tmp_path):
    exe, path = profile(tmp_path)
    with pytest.raises(ValueError, match='scope'):
        calibrated_estimate(path, exe, 2000, 64, 1024)


def test_large_image_workspace_is_included():
    assert native_estimate(31, 3, 1024)['planning_gib'] > native_estimate(31, 3, 128)['planning_gib']


@pytest.mark.parametrize('n,l,p', [(0, 1, 16), (1, -1, 16), (1, 1, -1), (True, 1, 16)])
def test_invalid_native_dimensions_rejected(n, l, p):
    with pytest.raises(ValueError):
        native_estimate(n, l, p)
