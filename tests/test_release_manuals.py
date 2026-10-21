from pathlib import Path

import pytest

from far_field_standalone.build_release import copy_release_manuals


def populate(root):
    for name in ('far_field_standalone/manual_assets/controls.png',
                 'docs/showcase_20261017/assets/far.png',
                 'docs/showcase_20261017/case_manifest.json',
                 'output/pdf/WaferScope_中文操作手册.pdf',
                 'output/pdf/WaferScope_功能与案例手册.pdf'):
        path = root/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b'fixture')


def test_manuals_preserve_assets_and_pdf_bytes(tmp_path):
    root, app = tmp_path/'source', tmp_path/'bundle'
    populate(root)
    copy_release_manuals(app, root)
    for name in ('docs/far_field_standalone/manual_assets/controls.png',
                 'docs/showcase_20261017/assets/far.png',
                 'docs/showcase_20261017/case_manifest.json',
                 'docs/WaferScope_中文操作手册.pdf',
                 'docs/WaferScope_功能与案例手册.pdf'):
        assert (app/name).read_bytes() == b'fixture'


def test_missing_manual_prevents_incomplete_release(tmp_path):
    populate(tmp_path)
    (tmp_path/'output/pdf/WaferScope_中文操作手册.pdf').unlink()
    with pytest.raises(FileNotFoundError, match='Missing release manual'):
        copy_release_manuals(tmp_path/'bundle', tmp_path)
