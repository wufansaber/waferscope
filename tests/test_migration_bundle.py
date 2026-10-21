import json
import zipfile

import pytest

from tools.build_migration_bundle import PREFIX, allowed, verify_bundle, write_bundle


@pytest.mark.parametrize('name', ['build/env/python.exe', '.vscode/settings.json', '.env',
    'docs/../private.txt', '/outside.txt', 'keys/private.pem', 'x\\outside.txt'])
def test_package_excludes_local_or_unsafe_paths(name):
    assert not allowed(name)


def test_package_is_hashed_and_omissions_explicit(tmp_path):
    source, target = tmp_path/'source.zip', tmp_path/'migration.zip'
    with zipfile.ZipFile(source, 'w') as archive:
        archive.writestr('docs/PROJECT_MIGRATION_CN.md', '# Migration\n')
        archive.writestr('rcwa_prototype/requirements.txt', 'numpy==2.4.6\n')
        archive.writestr('.vscode/settings.json', 'personal')
        archive.writestr('build/cache.npz', 'cache')
    manifest = write_bundle(source, target, dict(source_commit='test', source_tree='tree'),
        dict(rcwa=dict(python='3.11')), [dict(path='output/a/field.ngs', bytes=100)])
    assert len(manifest['files']) == 3
    assert set(manifest['excluded_tracked_paths']) == {'.vscode/settings.json', 'build/cache.npz'}
    assert manifest['omitted_local_outputs'][0]['path'] == 'output/a/field.ngs'
    assert verify_bundle(target) == manifest
    with zipfile.ZipFile(target) as archive:
        assert json.loads(archive.read(PREFIX+'/ENVIRONMENT_SNAPSHOTS.json'))['rcwa']['python'] == '3.11'
    with pytest.raises(FileExistsError):
        write_bundle(source, target, {}, {}, [])


def test_package_rejects_unlisted_payload(tmp_path):
    source, target = tmp_path/'source.zip', tmp_path/'migration.zip'
    with zipfile.ZipFile(source, 'w') as archive:
        archive.writestr('main.py', 'print(1)')
    write_bundle(source, target, {}, {}, [])
    with zipfile.ZipFile(target, 'a') as archive:
        archive.writestr(PREFIX+'/unexpected.txt', 'not listed')
    with pytest.raises(ValueError, match='Unlisted'):
        verify_bundle(target)
