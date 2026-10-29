"""Summarize saved S4 benchmarks without treating coarse particle trials as certified."""
import argparse
import json
from pathlib import Path

import numpy as np

from rcwa_prototype.comparison import compare_pair, field_change
from rcwa_prototype.geometry import Scene
from s4_prototype.validation import nist_comparison


def load(folder):
    info = json.loads((folder/'report.json').read_text(encoding='utf-8'))
    with np.load(folder/'modes.npz', allow_pickle=False) as saved:
        return info, {k: saved[k] for k in saved.files}


def summarize(root, oracle):
    trials, benchmarks, pairs, resources = {}, [], [], []
    for path in sorted(root.glob('*/report.json')):
        info, data = load(path.parent)
        if info.get('backend') != 'S4':
            continue
        name = path.parent.name
        trials[name] = (info, data)
        scene = Scene(**info['scene'])
        if not info['particle'] and scene.kind == 'grating' and info['polarization'] == 'p':
            benchmarks.append(dict(trial=name, one_dimensional=info['one_dimensional'],
                factorization=info['factorization'], actual_modes=info['actual_modes'],
                stripe_rule=info.get('stripe_rule', False),
                epsilon_method=info.get('epsilon_method', 'analytic'),
                retained_order_extent=info.get('retained_order_extent'),
                **nist_comparison(scene, data, oracle)))
        usage = path.parent.with_suffix('.resources.json')
        if usage.exists():
            saved = json.loads(usage.read_text(encoding='utf-8'))
            resources.append(dict(trial=name, **{key: saved[key] for key in (
                'status', 'elapsed_seconds', 'peak_accounted_gib', 'minimum_available_gib')}))
    for name, (info, sample) in trials.items():
        reference_name = name.replace('particle_', 'reference_', 1)
        if not name.startswith('particle_') or reference_name not in trials:
            continue
        reference_info, reference = trials[reference_name]
        if any(reference_info[key] != info[key] for key in (
                'scene', 'polarization', 'one_dimensional', 'factorization', 'resolution', 'executable_sha256')):
            raise ValueError('Incompatible sample/reference pair: '+name)
        for key, default in (('stripe_rule', False), ('epsilon_method', 'analytic')):
            if reference_info.get(key, default) != info.get(key, default):
                raise ValueError('Incompatible sample/reference material rule: '+name)
        scene = Scene(**info['scene'])
        pairs.append(dict(trial=name, actual_modes=info['actual_modes'], factorization=info['factorization'],
                          stripe_rule=info.get('stripe_rule', False),
                          harmonics=[scene.harmonics_x, scene.harmonics_y],
                          **compare_pair(reference, sample, scene)))
    changes = []
    for factorization in ('default', 'normal'):
        group = sorted([p for p in pairs if p['factorization'] == factorization], key=lambda p: p['actual_modes'])
        for coarse, fine in zip(group, group[1:]):
            delta = []
            for trial in (coarse, fine):
                name = trial['trial']
                _, sample = trials[name]
                _, reference = trials[name.replace('particle_', 'reference_', 1)]
                delta.append(dict(sample, electric_reflected_z0=
                    sample['electric_reflected_z0']-reference['electric_reflected_z0']))
            scene = Scene(**trials[fine['trial']][0]['scene'])
            old_scene = trials[coarse['trial']][0]['scene']
            changed = [key for key, value in vars(scene).items() if value != old_scene[key]]
            channel_changes = []
            for a, b in zip(coarse['channels'], fine['channels']):
                old, new = a['difference_field_efficiency'], b['difference_field_efficiency']
                channel_changes.append(old/new-1 if new else None)
            changes.append(dict(coarse=coarse['trial'], fine=fine['trial'],
                                changed_scene_parameters=changed,
                                channel_difference_field_efficiency_relative_changes=channel_changes,
                                **field_change(*delta, scene.incidence_deg)))
    for pair in pairs:
        a, b = [ch['difference_field_efficiency'] for ch in pair['channels'][:2]]
        pair['ch12_relative_asymmetry'] = abs(a-b)/((a+b)/2) if a+b else None
    passing = [b for b in benchmarks if
        abs(b['relative_total_reflectance_error']) < .01 and
        b['power_weighted_relative_complex_field_error'] < .02 and
        all(abs(row['relative_power_error']) < .05 for row in b['per_order'])]
    result = dict(status=('BARE_GRATING_VALIDATED_PARTICLE_NOT_CONVERGED' if passing else
                          'VALIDATION_INCOMPLETE'),
        benchmarks=benchmarks, pairs=pairs, truncation_changes=changes, resources=resources,
        passing_background_trials=[b['trial'] for b in passing],
        limitations=['One-dimensional bare-grating validation does not certify a 3D particle.',
            '3D particle truncation, sphere slicing and supercell separation are not converged.',
            'Legacy normal-vector factorization without the stripe adapter remains diagnostic only.',
            'The stripe adapter is a documented project modification, not an unmodified official S4 feature.',
            'Periodic array is not an isolated particle or a finite grating.',
            'Illustrative material constants are not calibrated at 266 nm.',
            'No camera image or quantitative SNR is certified by this report.'])
    (root/'summary.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    write_markdown(root, result)
    return result


def write_markdown(root, result):
    rows = ['# S4 validation results', '', 'Status: '+result['status'], '',
        '## Bare Grating vs NIST', '',
        '| Trial | Modes | Total R relative error | Weighted complex E error |',
        '|---|---:|---:|---:|']
    for b in result['benchmarks']:
        rows.append('| {} | {} | {:.6g}% | {:.6g}% |'.format(b['trial'], b['actual_modes'],
            b['relative_total_reflectance_error']*100, b['power_weighted_relative_complex_field_error']*100))
    rows += ['', '## Particle Trials (Not Converged)', '',
        '| Trial | CH | Sample efficiency | Reference efficiency | Difference-field efficiency | Actual power difference |',
        '|---|---:|---:|---:|---:|---:|']
    for pair in result['pairs']:
        for i, ch in enumerate(pair['channels'], 1):
            rows.append('| {} | {} | {:.6g} | {:.6g} | {:.6g} | {:.6g} |'.format(pair['trial'], i,
                ch['sample_efficiency'], ch['reference_efficiency'],
                ch['difference_field_efficiency'], ch['actual_power_difference']))
    rows += ['', 'All three channels: NA=0.6; tilts=(+45,-45,0) deg. Efficiencies are fractions of',
        'plane-wave power incident on the periodic supercell. Difference-field intensity is NOT sample intensity minus reference intensity.',
        '', '## Truncation Changes', '']
    rows += ['- {} -> {}: {:.6g}% weighted complex difference-field change.'.format(
        c['coarse'], c['fine'], c.get('relative_power_weighted_complex_field_change', float('nan'))*100)
        for c in result['truncation_changes']]
    rows += ['', '## Limitations', '']+['- '+item for item in result['limitations']]
    (root/'VALIDATION.md').write_text('\n'.join(rows)+'\n', encoding='utf-8')


def plot(root, result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    candidates = [b for b in result['benchmarks'] if b['factorization'] == 'default']
    if not candidates:
        return
    best = min(candidates, key=lambda b: b['power_weighted_relative_complex_field_error'])
    rows = best['per_order']
    x = np.arange(len(rows))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout='constrained')
    axes[0].bar(x-.17, [r['nist_efficiency'] for r in rows], .34, label='NIST RCW', color='#007D8A')
    axes[0].bar(x+.17, [r['s4_efficiency'] for r in rows], .34, label='S4', color='#CE6939')
    axes[0].set(yscale='log', xticks=x, xticklabels=[r['order'] for r in rows],
        xlabel='Diffraction order m', ylabel='Reflected power / incident power',
        title='Bare grating: NIST vs S4 ({} modes)'.format(best['actual_modes']))
    axes[0].legend(frameon=False)
    variants = sorted({(b['factorization'], b.get('stripe_rule', False), b.get('epsilon_method', 'analytic'))
        for b in result['benchmarks'] if not b['one_dimensional'] and b['factorization'] in ('default', 'normal')})
    colors = ['#007D8A', '#CE6939', '#555555', '#337A45']
    for index, (factorization, stripe, method) in enumerate(variants):
        group = sorted([b for b in result['benchmarks'] if not b['one_dimensional'] and
            (b['factorization'], b.get('stripe_rule', False), b.get('epsilon_method', 'analytic')) ==
            (factorization, stripe, method)], key=lambda b: b['actual_modes'])
        label = ('S4 + stripe adapter' if stripe else 'S4 legacy')+' / '+factorization+' / '+method
        axes[1].plot([b['actual_modes'] for b in group],
            [(1+b['relative_total_reflectance_error'])*best['oracle_inputs']['reflectance'] for b in group],
            'o-', label=label, color=colors[index % len(colors)])
    axes[1].axhline(best['oracle_inputs']['reflectance'], color='#333333', linestyle='--', label='NIST reference')
    axes[1].set(xlabel='Total modes in the 2D supercell basis', ylabel='Total reflectance',
                title='2D basis diagnostics (different cutoffs)')
    axes[1].legend(frameon=False)
    for axis in axes:
        axis.grid(axis='y', alpha=.2)
        axis.set_axisbelow(True)
    fig.savefig(root/'nist_s4_validation.png', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('output/s4_validation_20261022'))
    parser.add_argument('--nist', type=Path, default=Path('output/rcwa_grating_20261022/nist365.npz'))
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    result = summarize(args.root, args.nist)
    if args.plot:
        plot(args.root, result)
    print(json.dumps(dict(status=result['status'], comparisons=len(result['benchmarks']),
        truncation_changes=result['truncation_changes']), indent=2))


if __name__ == '__main__':
    main()
