import hashlib
import json
import zipfile

import pytest

from tools.restore_s4_dependencies import restore


@pytest.mark.parametrize('matches', [True, False])
def test_wheel_dll_is_restored_only_when_it_matches_build(tmp_path, matches):
    target = tmp_path/'native'
    target.mkdir()
    original = b'pinned test DLL'
    (target/'build_manifest.json').write_text(json.dumps(dict(
        blas_sha256=hashlib.sha256(original).hexdigest())), encoding='utf-8')
    wheel = tmp_path/'openblas.whl'
    with zipfile.ZipFile(wheel, 'w') as archive:
        archive.writestr('ngsolve_openblas/openblas.dll', original if matches else b'wrong DLL')
    if matches:
        restore(wheel, target)
        assert (target/'openblas.dll').read_bytes() == original
    else:
        with pytest.raises(ValueError, match='does not match'):
            restore(wheel, target)
        assert not (target/'openblas.dll').exists()
