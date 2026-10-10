"""Run one FEM job with sampled memory/time limits and child-process cleanup."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time


def stop_reason(memory_gib, available_gib, elapsed, maximum, reserve, timeout):
    if memory_gib > maximum:
        return 'memory_limit'
    if available_gib < reserve:
        return 'system_memory_reserve'
    if elapsed > timeout:
        return 'timeout'
    return None


def _wait_for_owned(processes, timeout=10):
    """Wait only for living processes; Windows may deny waits on exited PIDs."""
    import psutil

    def still_running(candidates):
        result = []
        for process in candidates:
            try:
                if process.is_running():
                    result.append(process)
            except psutil.NoSuchProcess:
                pass
            except psutil.AccessDenied:
                result.append(process)  # Unknown is not silently treated as gone.
        return result

    candidates = still_running(processes)
    if not candidates:
        return []
    try:
        _, alive = psutil.wait_procs(candidates, timeout=timeout)
        return [process.pid for process in still_running(alive)]
    except psutil.AccessDenied:
        time.sleep(min(timeout, 1))
        return [process.pid for process in still_running(candidates)]


def run_guarded(command, report_path, maximum=10., reserve=3., timeout=1800., cancel_file=None):
    import psutil
    for value in (maximum, reserve, timeout):
        if not math.isfinite(value) or value <= 0:
            raise ValueError('Memory, reserve and timeout limits must be positive and finite')
    path = Path(report_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if psutil.virtual_memory().available/2**30 < reserve+.25:
        raise RuntimeError('Insufficient available memory to start the guarded job')
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1',
               PYTHONUNBUFFERED='1')
    start = time.monotonic()
    child = subprocess.Popen(command, env=env)
    root = psutil.Process(child.pid)
    owned = {child.pid: root}
    report = dict(command=command, limits=dict(maximum_gib=maximum, reserve_gib=reserve,
                  timeout_seconds=timeout), peak_rss_gib=0., peak_accounted_gib=0.,
                  minimum_available_gib=psutil.virtual_memory().available/2**30,
                  status='running', sampling_seconds=.25,
                  note='Sampled limits, not an OS hard memory cap. Windows accounting includes private committed bytes; shared RSS may be counted twice.')
    try:
        while child.poll() is None:
            try:
                for process in root.children(recursive=True):
                    owned[process.pid] = process
            except psutil.NoSuchProcess:
                pass
            rss, accounted = 0, 0
            for process in owned.values():
                try:
                    info = process.memory_info()
                    rss += info.rss
                    accounted += max(info.rss, getattr(info, 'private', info.rss))
                except psutil.NoSuchProcess:
                    pass
            available = psutil.virtual_memory().available/2**30
            elapsed = time.monotonic()-start
            report['peak_rss_gib'] = max(report['peak_rss_gib'], rss/2**30)
            report['peak_accounted_gib'] = max(report['peak_accounted_gib'], accounted/2**30)
            report['minimum_available_gib'] = min(report['minimum_available_gib'], available)
            reason = ('cancelled' if cancel_file and Path(cancel_file).exists() else
                      stop_reason(accounted/2**30, available, elapsed, maximum, reserve, timeout))
            if reason:
                report['status'] = 'stopped_'+reason
                print('Resource guard stopping its own job: '+reason, flush=True)
                break
            time.sleep(.25)
        if report['status'] == 'running':
            report['status'] = 'completed' if child.returncode == 0 else 'failed'
    except BaseException:
        report['status'] = 'interrupted_or_monitor_error'
        raise
    finally:
        # Process objects check creation identity; never terminate unrelated named processes.
        for process in reversed(list(owned.values())):
            try:
                if process.is_running():
                    process.kill()
            except psutil.NoSuchProcess:
                pass
        alive = _wait_for_owned(list(owned.values()), timeout=10)
        report.update(elapsed_seconds=time.monotonic()-start, returncode=child.wait(timeout=10),
                      surviving_owned_pids=alive)
        if alive:
            report['status'] = 'cleanup_failed'
        path.write_text(json.dumps(report, indent=2), encoding='utf-8')
        print(json.dumps(report, indent=2), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', required=True)
    parser.add_argument('--max-memory-gib', type=float, default=10)
    parser.add_argument('--reserve-gib', type=float, default=3)
    parser.add_argument('--timeout-seconds', type=float, default=1800)
    parser.add_argument('--cancel-file')
    parser.add_argument('module', choices=('physical_audit', 'imaging_run', 'rough_run', 'theory_audit',
                                           'near_audit', 'component_compare', 'haze_ensemble',
                                           'haze_convergence', 'whole_haze', 'hybrid_run',
                                           'local_spm_validation', 'hybrid_range_validation', 'surface_mesh_audit',
                                           'illumination_definition_audit', 'range_refinement'))
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    result = run_guarded([sys.executable, '-m', 'fullwave_prototype.'+args.module]+args.arguments,
                         args.report, args.max_memory_gib, args.reserve_gib, args.timeout_seconds,
                         args.cancel_file)
    raise SystemExit(0 if result['status'] == 'completed' else 1)


if __name__ == '__main__':
    main()
