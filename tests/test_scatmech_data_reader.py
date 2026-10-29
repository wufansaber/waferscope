"""Independent SCATMECH reader supports new exports and numeric legacy NPZ."""
import numpy as np
import pytest

from result_archive import save_result_archive
from far_field_standalone.read_scatmech_data import read_file


def test_reader_checks_sha_and_rebuilds_all_channels(tmp_path):
    path = tmp_path/'far_field_data.npz'
    channels = [dict(sampled_dsc_um2_sr=np.full((2, 3), i*1e-9),
                     aperture_mask=np.ones((2, 3), bool)) for i in (1, 2, 3)]
    save_result_archive(path, dict(plot_kind='objective_collection', channel_pupils=channels),
                        source=dict(jones_dsc_um=np.ones((2, 3, 2, 2), complex)))
    result, source, summary = read_file(path)
    assert summary['manifest_verified'] is True
    assert summary['plot_kind'] == 'objective_collection'
    assert len(result['channel_pupils']) == 3
    np.testing.assert_array_equal(result['channel_pupils'][2]['sampled_dsc_um2_sr'], channels[2]['sampled_dsc_um2_sr'])
    assert np.iscomplexobj(source['jones_dsc_um'])
    with path.open('ab') as stream:
        stream.write(b'changed')
    with pytest.raises(ValueError, match='SHA256'):
        read_file(path)


def test_reader_legacy_and_no_manifest_do_not_claim_checksum_verification(tmp_path):
    path = tmp_path/'old.npz'
    original = np.array([[1+2j, 3-4j]])
    np.savez_compressed(path, field=original)
    result, source, summary = read_file(path)
    assert source is None
    assert summary['manifest_verified'] is None
    assert summary['format'] == 'legacy-flat-npz'
    np.testing.assert_array_equal(result['field'], original)
    assert summary['arrays']['field'] == dict(shape=[1, 2], dtype=original.dtype.str, complex=True)


def test_reader_modern_no_manifest_still_reconstructs(tmp_path):
    path = tmp_path/'modern.npz'
    save_result_archive(path, dict(parameters=dict(wavelength_nm=266.), value=np.ones(3)))
    path.with_suffix('.manifest.json').unlink()
    result, _, summary = read_file(path)
    assert result['parameters']['wavelength_nm'] == 266.
    assert summary['manifest_verified'] is None
