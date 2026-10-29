"""Guarded periodic scene CLI used by the Qt page and headless workflows."""
import argparse
import json
from pathlib import Path
import sys

import psutil

from fullwave_prototype.guarded_run import run_guarded
from rcwa_prototype.geometry import matrix_memory_estimate
from .config import SceneConfig
from .pipeline import DEFAULT_EXECUTABLE, reimage_scene, run_scene


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--executable', type=Path, default=DEFAULT_EXECUTABLE)
    parser.add_argument('--max-memory-gib', type=float, default=8.)
    parser.add_argument('--timeout', type=float, default=3600.)
    parser.add_argument('--cancel-file', type=Path)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--reuse-source', type=Path)
    parser.add_argument('--memory-profile', type=Path)
    args = parser.parse_args()
    config = SceneConfig.from_dict(json.loads(args.config.read_text(encoding='utf-8')))
    if args.worker:
        if args.reuse_source:
            reimage_scene(config, args.reuse_source, args.output)
        else:
            run_scene(config, args.output, args.executable, args.timeout)
        return
    estimate = matrix_memory_estimate(config.scene.requested_modes, config.scene.slices+1)
    if args.reuse_source:
        estimate['planning_gib'] = .4+config.image_samples**2*16*40/2**30
    elif args.memory_profile:
        from .resources import calibrated_estimate
        estimate = calibrated_estimate(args.memory_profile, args.executable,
            config.scene.requested_modes, config.scene.slices+1, config.image_samples)
    available = psutil.virtual_memory().available/2**30
    if estimate['planning_gib'] > min(args.max_memory_gib, available-3.):
        parser.error('Planning {:.3f} GiB exceeds memory budget {:.3f} GiB (3 GiB system reserve)'.format(
            estimate['planning_gib'], max(0., min(args.max_memory_gib, available-3.))))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    from .pipeline import write_json
    write_json(args.output.with_suffix('.admission.json'), dict(planning=estimate,
        available_gib=available, maximum_gib=args.max_memory_gib, system_reserve_gib=3.))
    result = run_guarded([sys.executable, '-m', 'grating_scene.run', '--worker']+sys.argv[1:],
        args.output.with_suffix('.resources.json'), maximum=args.max_memory_gib, reserve=3.,
        timeout=args.timeout, cancel_file=args.cancel_file)
    if result['status'] != 'completed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
