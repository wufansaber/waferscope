"""Raw archives keep nested channels, phases, inputs and legacy numeric keys."""
from dataclasses import dataclass
import hashlib
import json
import numpy as np
import pytest
from result_archive import save_result_archive, load_result_archive, verify_result_archive
from fullwave_prototype.source_snapshot import snapshot_sources, HYBRID_FILES


@dataclass
class Inputs:
    wavelength_nm: float = 266.


def test_nested_raw_roundtrip(tmp_path):
    field = np.array([1+2j, -3j])
    raw = np.array([1e-15, np.nan])
    original = dict(parameters=Inputs(), channels=[dict(field=field, raw=raw)]*3,
                    count=np.int64(3), mask=np.array([True, False]),
                    angles=(0., np.pi), missing=float('nan'))
    path = tmp_path/'result.npz'
    manifest = save_result_archive(path, original, source=dict(raw=raw),
                                   legacy_arrays=dict(ch1_field=field))
    result, source = load_result_archive(path)
    assert result['parameters']['wavelength_nm'] == 266.
    assert len(result['channels']) == 3
    np.testing.assert_array_equal(result['channels'][2]['field'], field)
    np.testing.assert_array_equal(source['raw'], raw)
    assert isinstance(result['angles'], tuple)
    assert np.isnan(result['missing'])
    assert manifest['sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()
    with np.load(path, allow_pickle=False) as saved:
        np.testing.assert_array_equal(saved['ch1_field'], field)
        assert all(not saved[key].dtype.hasobject for key in saved.files)
    assert json.loads(path.with_suffix('.manifest.json').read_text()) == manifest


def test_reject_object_preserves_existing_result(tmp_path):
    path = tmp_path/'result.npz'
    save_result_archive(path, dict(raw=np.ones(3)))
    before = path.read_bytes()
    with pytest.raises(ValueError, match='Object arrays'):
        save_result_archive(path, dict(raw=np.array([{}], dtype=object)))
    assert path.read_bytes() == before


def test_detect_corrupt_inventory(tmp_path):
    path = tmp_path/'result.npz'
    save_result_archive(path, dict(raw=np.ones(3)))
    with np.load(path, allow_pickle=False) as saved:
        arrays = {key: saved[key] for key in saved.files}
    arrays['raw'] = np.zeros(4)
    np.savez_compressed(path, **arrays)
    with pytest.raises(ValueError, match='specification differs'):
        verify_result_archive(path)


@pytest.mark.parametrize('hybrid', [False, True])
def test_portable_source_snapshot(tmp_path, hybrid):
    source, result = tmp_path/'source', tmp_path/'result'
    for name in HYBRID_FILES if hybrid else ('report.json', 'paired_fields.npz'):
        path = source/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(name.encode('ascii'))
    snapshot_sources(source, result, hybrid=hybrid)
    snapshot_sources(source, result, hybrid=hybrid)
    manifest = json.loads((result/'source_manifest.json').read_text())
    for name, spec in manifest['files'].items():
        assert hashlib.sha256((result/name).read_bytes()).hexdigest() == spec['sha256']
    (result/name).write_bytes(b'different')
    with pytest.raises(ValueError, match='Refusing to overwrite'):
        snapshot_sources(source, result, hybrid=hybrid)
