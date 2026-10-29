"""Shared catalog for UI presets and complete, portable demonstration records."""
import json
from pathlib import Path

from .config import SceneConfig

EXAMPLES = (
    dict(key='contact_preview', title='100 nm sphere + grating | Preview', kind='scene',
         source='output/grating_scene_20261027/contact_preview'),
    dict(key='contact_near_far_image', title='100 nm sphere + grating | Near / Far / Image', kind='scene',
         source='output/grating_scene_20261027/native_controls/joint'),
    dict(key='contact_camera_40x', title='100 nm sphere + grating | 40x camera', kind='camera',
         source='output/grating_scene_20261027/camera_gemi_area'))


def example_record(key):
    for row in EXAMPLES:
        if row['key'] == key:
            return row
    raise ValueError('Unknown grating example: '+str(key))


def example_folder(root, key):
    row = example_record(key)
    root = Path(root)
    bundled = root/'examples/grating'/key
    return bundled if bundled.is_dir() else root/row['source']


def example_inputs(root, key):
    row = example_record(key)
    folder = example_folder(root, key)
    optical = folder/'optical_scene' if row['kind'] == 'camera' else folder
    config = SceneConfig.from_dict(json.loads((optical/'config.json').read_text(encoding='utf-8')))
    camera = json.loads((folder/'camera_settings.json').read_text(encoding='utf-8')) if row['kind'] == 'camera' else None
    return config, camera
