import pytest

from tools.publish_git_snapshot import commit_payload, identity


def test_git_identity_preserves_timezone():
    value = identity('Example User <example@example.invalid> 0 +0800')
    assert value['date'] == '1970-01-01T08:00:00+08:00'
    assert value['name'] == 'Example User'


def test_commit_metadata_preserves_message_and_single_parent():
    raw = (b'tree tree123\nparent parent123\nauthor Example <e@example.invalid> 0 +0000\n'
           b'committer Example <e@example.invalid> 0 +0000\n\nMigration snapshot\n')
    data = commit_payload(raw, 'tree123')
    assert data['message'] == 'Migration snapshot\n'
    assert data['parents'] == ['parent123']
    assert data['committer']['date'].endswith('+00:00')
    with pytest.raises(ValueError, match='identical'):
        commit_payload(raw, 'wrong_tree')


def test_signed_commit_is_not_silently_rewritten():
    raw = (b'tree t\nparent p\nauthor Example <e@example.invalid> 0 +0000\n'
           b'committer Example <e@example.invalid> 0 +0000\ngpgsig signature\n\nMessage\n')
    with pytest.raises(ValueError, match='unsigned'):
        commit_payload(raw, 't')
