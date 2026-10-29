"""Physical quadrature, absolute complex radiation and complete bundle contracts."""
from dataclasses import asdict
import hashlib
import json

import numpy as np
import pytest

from rcwa_prototype import angular_distribution as angular
from rcwa_prototype.geometry import Scene


def test_hemisphere_solid_angle_and_first_moment():
    arrays = angular.hemisphere_grid(24, 72)
    w, n = arrays['solid_angle_weights'], arrays['directions']
    assert w.sum() == pytest.approx(2*np.pi, rel=1e-14)
    np.testing.assert_allclose(np.linalg.norm(n, axis=-1), 1., atol=1e-14)
    np.testing.assert_allclose(np.sum(w[:, :, None]*n, axis=(0, 1)), [0, 0, np.pi], atol=2e-14)
    assert np.all(n[:, :, 2] > 0)
    assert np.all(np.diff(arrays['theta_deg']) > 0)


@pytest.mark.parametrize('polar,azimuthal', [(True, 8), (4., 8), (3, 8), (4, 7), (257, 8), (4, 1025)])
def test_bad_grid_rejected(polar, azimuthal):
    with pytest.raises(ValueError):
        angular.hemisphere_grid(polar, azimuthal)


def test_dipole_pattern_absolute_integral_and_phase(monkeypatch):
    def dipole(source, scene, ray, **kwargs):
        # Unit Z dipole has DSC=sin(theta)^2 and upper-hemisphere integral 4*pi/3.
        return (np.array([0., 0., 1.])-ray*ray[2])*np.exp(.37j)*2e-5
    monkeypatch.setattr(angular, 'project_ray', dipole)
    arrays = angular.evaluate({}, Scene(), 12, 36, 20)
    assert np.sum(arrays['solid_angle_weights']*arrays['dsc_um2_per_sr']) == pytest.approx(4e-10*4*np.pi/3, rel=1e-14)
    np.testing.assert_allclose(np.sum(arrays['directions']*arrays['far_electric_um'], axis=-1), 0, atol=1e-20)
    assert arrays['far_electric_um'][0, 0, 2].imag > 0


def fixtures(tmp_path):
    scene = Scene(grating_method='analytic_li')
    paths = []
    for name, particle in (('particle', True), ('bare', False)):
        folder = tmp_path/name
        folder.mkdir()
        np.savez_compressed(folder/'modes.npz', source_electric_fourier=np.zeros((1, 1, 3)),
            electric_reflected_z0=np.ones((1, 3)), directions=np.array([[0., 0., 1.]]), orders=np.array([[0, 0]]))
        info = dict(scene=asdict(scene), particle=particle, current_z_order=8,
            backend='grcwa', backend_version='0.1.3', coordinates='z up',
            raw_npz_sha256=hashlib.sha256((folder/'modes.npz').read_bytes()).hexdigest())
        (folder/'report.json').write_text(json.dumps(info), encoding='utf-8')
        paths.append(folder)
    return paths


def mock_gate(monkeypatch, error=0.):
    monkeypatch.setattr(angular, 'modal_reconstruction', lambda *a: (np.zeros((1, 3)),
        dict(comparable=True, relative_power_weighted_complex_field_change=error, relative_total_power_difference=0.)))
    monkeypatch.setattr(angular, 'project_ray', lambda *a, **kw: np.array([1e-5, 2e-5j, 0.]))


def test_bundle_roundtrip_export_data_and_integrity(tmp_path, monkeypatch):
    source, reference = fixtures(tmp_path)
    mock_gate(monkeypatch)
    output = tmp_path/'result'
    angular.run(source, reference, output, 4, 8, 2)
    report, arrays = angular.load_distribution(output)
    assert report['hemisphere_cross_section_um2'] == pytest.approx(5e-10*2*np.pi)
    assert report['isolated_defect_certified'] is False
    assert arrays['far_electric_um'].shape == (4, 8, 3)
    csv = np.loadtxt(output/'angular_distribution.csv', delimiter=',', skiprows=1)
    np.testing.assert_allclose(csv[:, -1], arrays['dsc_um2_per_sr'].ravel())
    np.testing.assert_allclose(csv[:, 6:9]+1j*csv[:, 9:12], arrays['far_electric_um'].reshape(-1, 3))
    with pytest.raises(FileExistsError):
        angular.run(source, reference, output, 4, 8, 2)
    (output/'angular_distribution.csv').write_text('corrupt', encoding='utf-8')
    with pytest.raises(ValueError, match='hash'):
        angular.load_distribution(output)


def test_failed_gate_has_no_completion_or_radiation(tmp_path, monkeypatch):
    source, reference = fixtures(tmp_path)
    mock_gate(monkeypatch, error=.01)
    with pytest.raises(RuntimeError, match='gate failed'):
        angular.run(source, reference, tmp_path/'result', 4, 8, 2)
    assert not (tmp_path/'result').exists()


def test_nonfinite_radiation_rejected(monkeypatch):
    monkeypatch.setattr(angular, 'project_ray', lambda *a, **kw: np.array([np.nan, 0, 0]))
    with pytest.raises(ValueError, match='Nonfinite'):
        angular.evaluate({}, Scene(), 4, 8, 2)


def test_mismatched_metadata_rejected(tmp_path, monkeypatch):
    source, reference = fixtures(tmp_path)
    info = json.loads((reference/'report.json').read_text(encoding='utf-8'))
    info['coordinates'] = 'wrong phase reference'
    (reference/'report.json').write_text(json.dumps(info), encoding='utf-8')
    with pytest.raises(ValueError, match='metadata'):
        angular.run(source, reference, tmp_path/'result', 4, 8, 2)
