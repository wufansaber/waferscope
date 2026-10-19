import hashlib
import json

import numpy as np
import pytest

from tools.build_showcase_cases import load_dense_far


def make_cache(tmp_path):
    source = tmp_path/'source'; assets = tmp_path/'assets'
    source.mkdir(); assets.mkdir()
    metadata = source/'hybrid_source_report.json'
    metadata.write_text('{"physical_source": 1}', encoding='utf-8')
    report = dict(source_metadata_sha256=hashlib.sha256(metadata.read_bytes()).hexdigest())
    (assets/'case01_dense_far_report.json').write_text(json.dumps(report), encoding='utf-8')
    axis = np.linspace(-1.,1.,401)
    x,y = np.meshgrid(axis,axis); mask = x*x+y*y<1-1e-12
    values = np.where(mask, 1., np.nan)
    np.savez(assets/'case01_dense_far.npz',axis=axis,dsc_total=values,mask=mask)
    return source,assets,axis,mask,values


def test_dense_cache_uses_original_source_and_absolute_values(tmp_path):
    source,assets,axis,mask,values = make_cache(tmp_path)
    actual_axis,actual_values,_ = load_dense_far(assets,source)
    np.testing.assert_array_equal(actual_axis,axis)
    np.testing.assert_array_equal(actual_values,values)


def test_dense_cache_rejects_changed_physical_source(tmp_path):
    source,assets,*_ = make_cache(tmp_path)
    (source/'hybrid_source_report.json').write_text('{"physical_source": 2}',encoding='utf-8')
    with pytest.raises(ValueError,match='different physical source'):
        load_dense_far(assets,source)


@pytest.mark.parametrize('bad', ['negative','nonfinite','outside'])
def test_dense_cache_rejects_invalid_intensity(tmp_path,bad):
    source,assets,axis,mask,values = make_cache(tmp_path)
    if bad == 'outside':
        values[0,0]=1.
    else:
        values[200,200] = -1 if bad=='negative' else np.nan
    np.savez(assets/'case01_dense_far.npz',axis=axis,dsc_total=values,mask=mask)
    with pytest.raises(ValueError,match='Invalid dense'):
        load_dense_far(assets,source)
