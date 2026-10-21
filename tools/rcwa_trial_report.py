"""Rebuild the experimental convergence/resource report from guarded trials."""
import argparse
import json
from pathlib import Path

import numpy as np
import psutil

from rcwa_prototype.comparison import compare_pair, field_change
from rcwa_prototype.geometry import Scene


def read_trial(folder):
    info = json.loads((folder/'report.json').read_text(encoding='utf-8'))
    with np.load(folder/'modes.npz') as data:
        arrays = {k: data[k] for k in data.files}
    return info, arrays


def nist_comparison(root, trials):
    path = root/'nist245.npz'
    if not path.exists():
        return []
    with np.load(path) as oracle:
        expected = dict(zip(oracle['orders'].tolist(), oracle['efficiency'].tolist()))
    result = []
    for folder, info, data in trials:
        if info['geometry']['particle_present'] or info['scene']['kind'] != 'grating':
            continue
        p = info['scene']['periods_y']
        actual = {int(g[1]//p): float(v) for g, v in zip(data['orders'], data['reflectance_by_order'])
                  if g[0] == 0 and g[1] % p == 0}
        per_order = [dict(order=m, nist=e, rcwa=actual.get(m, 0.),
            relative_difference=(actual.get(m, 0.)/e-1)) for m, e in expected.items() if e > 1e-12]
        result.append(dict(trial=folder.name, actual_modes=info['actual_modes'],
            relative_total_reflectance_difference=info['reflectance']/sum(expected.values())-1,
            per_order=per_order))
    return result


def collect_pairs(trials):
    lookup = {folder.name: (info, data) for folder, info, data in trials}
    pairs = []
    for folder, info, sample in trials:
        if not folder.name.startswith('particle_'):
            continue
        name = folder.name.replace('particle_', 'reference_', 1)
        if name not in lookup:
            continue
        reference_info, reference = lookup[name]
        if reference_info['scene'] != info['scene']:
            raise ValueError('Mismatched reference scene: '+folder.name)
        pairs.append((folder.name, reference, sample, Scene(**info['scene'])))
    return sorted(pairs, key=lambda pair: pair[-1].modes)


def summarize(root):
    trials = [(folder, *read_trial(folder)) for folder in
              sorted(p.parent for p in root.glob('*/report.json'))]
    resources = []
    for folder, info, _ in trials:
        usage = json.loads(folder.with_suffix('.resources.json').read_text(encoding='utf-8'))
        resources.append(dict(trial=folder.name, actual_modes=info['actual_modes'],
            patterned_layers=info['patterned_layers'], reflectance=info['reflectance'],
            **{k: usage[k] for k in ('status', 'elapsed_seconds', 'peak_rss_gib',
                                    'peak_accounted_gib', 'minimum_available_gib')}))
    pairs = collect_pairs(trials)
    comparisons = [dict(trial=name, **compare_pair(ref, sample, scene))
                   for name, ref, sample, scene in pairs]
    changes = []
    for a, b in zip(pairs, pairs[1:]):
        if {k: v for k, v in vars(a[-1]).items() if k != 'modes'} != {
                k: v for k, v in vars(b[-1]).items() if k != 'modes'}:
            continue
        da = dict(a[2], electric_reflected_z0=a[2]['electric_reflected_z0']-a[1]['electric_reflected_z0'])
        db = dict(b[2], electric_reflected_z0=b[2]['electric_reflected_z0']-b[1]['electric_reflected_z0'])
        changes.append(dict(coarse=a[0], fine=b[0],
            difference_field=field_change(da, db, b[-1].incidence_deg)))
    memory = psutil.virtual_memory()
    particle_usage = [r for r in resources if r['trial'].startswith('particle_')]
    scaling = {}
    if particle_usage:
        largest = max(particle_usage, key=lambda row: row['actual_modes'])
        scaling = dict(based_on_trial=largest['trial'], patterned_layers=largest['patterned_layers'],
            estimates=[dict(total_modes=m, estimated_process_gib=largest['peak_accounted_gib']*
                (m/largest['actual_modes'])**2) for m in (1200, 2400, 4000)],
            note='Unmeasured quadratic extrapolation at fixed layers/backend; add OS/workspace headroom.')
    report = dict(status='EXPERIMENTAL_NOT_CONVERGED', resources=resources,
        hardware=dict(total_memory_gib=memory.total/2**30, available_at_report_gib=memory.available/2**30),
        memory_scaling=scaling,
        nist_comparisons=nist_comparison(root, trials), pair_comparisons=comparisons,
        truncation_changes=changes,
        limitations=['Periodic particle copies; isolated-defect limit not demonstrated.',
            'Backend is grcwa; S4 and AFMM have not been run.',
            'Legacy Si optical constants are not a calibrated 266 nm material input.',
            'Resource monitoring is sampled, not an OS hard allocation limit.',
            'Sparse order amplitudes have not been converted into a camera prediction.'])
    (root/'summary.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(summarize(args.root), indent=2))


if __name__ == '__main__':
    main()
