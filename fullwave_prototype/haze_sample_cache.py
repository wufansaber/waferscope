"""Accept only matching, successfully guarded haze samples for resumption."""
import json


def completed_sample(folder, request):
    required = ('input.json', 'resources.json', 'sample_report.json', 'haze_sources.npz')
    if not all((folder/name).is_file() for name in required):
        return False
    try:
        previous = json.loads((folder/'input.json').read_text(encoding='utf-8'))
        resources = json.loads((folder/'resources.json').read_text(encoding='utf-8'))
        return (previous == request and resources['status'] == 'completed'
                and resources['returncode'] == 0 and resources['surviving_owned_pids'] == [])
    except (OSError, ValueError, KeyError, TypeError):
        return False
