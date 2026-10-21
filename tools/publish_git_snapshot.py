"""Publish one exact local commit via authenticated gh API when Git HTTPS is blocked.

Creates a new codex branch only. Never rewrites an existing remote branch or reads tokens.
"""
import argparse
import base64
from datetime import datetime, timedelta, timezone
import json
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def git(*args):
    return subprocess.check_output(['git', '-c', 'safe.directory='+ROOT.as_posix(),
        '-C', str(ROOT), *args])


def identity(value):
    match = re.fullmatch(r'(.*) <([^<>]+)> (\d+) ([+-]\d{4})', value)
    if not match:
        raise ValueError('Unsupported Git identity')
    name, email, stamp, zone = match.groups()
    minutes = (int(zone[1:3])*60+int(zone[3:]))*(1 if zone[0] == '+' else -1)
    date = datetime.fromtimestamp(int(stamp), timezone(timedelta(minutes=minutes))).isoformat()
    return dict(name=name, email=email, date=date)


def commit_payload(raw, tree):
    headers, message = raw.decode('utf-8').split('\n\n', 1)
    pairs = [line.split(' ', 1) for line in headers.splitlines()]
    if any(key not in ('tree', 'parent', 'author', 'committer') for key, _ in pairs):
        raise ValueError('Only unsigned standard UTF-8 commits supported')
    fields = dict(pairs)
    parents = [value for key, value in pairs if key == 'parent']
    if len(parents) != 1 or fields['tree'] != tree:
        raise ValueError('Exactly one parent and an identical tree required')
    return dict(message=message, tree=tree, parents=parents,
        author=identity(fields['author']), committer=identity(fields['committer']))


def api(repo, endpoint, payload=None):
    command = ['gh', 'api', 'repos/'+repo+'/'+endpoint]
    if payload is None:
        return json.loads(subprocess.check_output(command, text=True, encoding='utf-8', timeout=180))
    with tempfile.TemporaryDirectory(prefix='github_api_', dir=ROOT/'build') as temporary:
        request = Path(temporary)/'request.json'
        request.write_text(json.dumps(payload), encoding='utf-8')
        return json.loads(subprocess.check_output(command+['-X', 'POST', '--input', str(request)],
            text=True, encoding='utf-8', timeout=180))


def publish(repo, branch):
    if repo != 'wufansaber/waferscope' or not branch.startswith('codex/'):
        raise ValueError('This migration uploader is scoped to wufansaber/waferscope and a new codex branch')
    git('check-ref-format', '--branch', branch)
    commit = git('rev-parse', 'HEAD').decode().strip()
    tree = git('rev-parse', 'HEAD^{tree}').decode().strip()
    payload = commit_payload(git('cat-file', 'commit', commit), tree)
    parent = payload['parents'][0]
    base = api(repo, 'git/commits/'+parent)
    existing = api(repo, 'git/matching-refs/heads/'+branch)
    if existing:
        if len(existing) == 1 and existing[0]['ref'] == 'refs/heads/'+branch and existing[0]['object']['sha'] == commit:
            return dict(repository=repo, branch=branch, commit=commit, status='already_published')
        raise ValueError('Remote branch already exists; refusing to overwrite it')
    blobs = {}
    for record in git('ls-tree', '-r', '-z', '--full-tree', 'HEAD').split(b'\0'):
        if record:
            header, path = record.split(b'\t', 1)
            mode, kind, sha = header.decode().split()
            blobs[path.decode('utf-8')] = (mode, kind, sha)
    changed = git('diff-tree', '-r', '--no-commit-id', '--no-renames', '--name-only', '-z', parent, commit)
    entries = []
    for raw in changed.split(b'\0'):
        if not raw:
            continue
        path = raw.decode('utf-8')
        if path not in blobs:
            entries.append(dict(path=path, mode='100644', type='blob', sha=None))
            continue
        mode, kind, sha = blobs[path]
        if kind != 'blob' or mode not in ('100644', '100755'):
            raise ValueError('Unsupported changed entry: '+path)
        data = git('cat-file', 'blob', sha)
        if len(data) > 50*1024**2:
            raise ValueError('Do not upload oversized cache: '+path)
        entry = dict(path=path, mode=mode, type='blob')
        try:
            content = data.decode('utf-8')
            if '\0' in content:
                raise UnicodeError('Binary data')
            entry['content'] = content
        except UnicodeError:
            response = api(repo, 'git/blobs', dict(content=base64.b64encode(data).decode('ascii'), encoding='base64'))
            if response['sha'] != sha:
                raise ValueError('GitHub blob hash differs: '+path)
            entry['sha'] = sha
        entries.append(entry)
        print('Prepared '+path, flush=True)
    created = api(repo, 'git/trees', dict(base_tree=base['tree']['sha'], tree=entries))
    if created['sha'] != tree:
        raise ValueError('Remote tree differs from local committed source; branch not created')
    created = api(repo, 'git/commits', payload)
    if created['sha'] != commit:
        raise ValueError('Remote commit differs from exact local commit; branch not created: '+created['sha'])
    api(repo, 'git/refs', dict(ref='refs/heads/'+branch, sha=commit))
    verified = api(repo, 'git/ref/heads/'+branch)
    if verified['object']['sha'] != commit:
        raise ValueError('Remote branch verification failed')
    return dict(repository=repo, branch=branch, commit=commit, source_tree=tree, status='published',
        url='https://github.com/'+repo+'/tree/'+branch)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--branch', required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    result = publish('wufansaber/waferscope', args.branch)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
