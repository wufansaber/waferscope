"""Audit explicit one-variable S4 convergence steps from retained complex fields."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from rcwa_prototype.comparison import compare_pair, field_change
from rcwa_prototype.geometry import Scene
from s4_prototype.solver import ROOT, sha256
from s4_prototype.validation import nist_comparison
from tools.s4_validation_report import load


SETTINGS = ('backend', 'polarization', 'one_dimensional', 'factorization',
            'resolution', 'executable_sha256')
NUMERICAL = ('harmonics_x', 'harmonics_y', 'slices')


def same_settings(a, b):
    for key in SETTINGS:
        if a[key] != b[key]:
            raise ValueError('Incompatible solver setting: '+key)
    for key, default in (('stripe_rule', False), ('epsilon_method', 'analytic')):
        if a.get(key, default) != b.get(key, default):
            raise ValueError('Incompatible material rule: '+key)


def load_checked(folder):
    info, data = load(folder)
    if sha256(folder/'modes.npz') != info['raw_npz_sha256']:
        raise ValueError('Raw NPZ hash mismatch: '+str(folder))
    if sha256(folder/'input.lua') != info['input_sha256']:
        raise ValueError('Lua input hash mismatch: '+str(folder))
    if not np.all(np.isfinite(data['electric_reflected_z0'])):
        raise ValueError('Nonfinite field: '+str(folder))
    return info, data


def prepare_pair(reference_info, reference, info, sample):
    same_settings(reference_info, info)
    if info['backend'] != 'S4' or not info['particle'] or reference_info['particle']:
        raise ValueError('Expected a bare-reference/S4 particle pair')
    # Bare cap is one homogeneous air layer, independent of sphere slice count.
    different = [k for k in info['scene'] if info['scene'][k] != reference_info['scene'][k]]
    if any(k != 'slices' for k in different):
        raise ValueError('Incompatible reference geometry: '+str(different))
    scene = Scene(**info['scene'])
    summary = compare_pair(reference, sample, scene)
    delta = dict(sample, electric_reflected_z0=
                 sample['electric_reflected_z0']-reference['electric_reflected_z0'])
    a, b = [ch['difference_field_efficiency'] for ch in summary['channels'][:2]]
    summary.update(actual_modes=info['actual_modes'], scene=asdict(scene),
        reflectance=info['reflectance'], reference_reflectance=reference_info['reflectance'],
        ch12_relative_asymmetry=abs(a-b)/((a+b)/2) if a+b else None,
        reference_reused_across_slice_count=bool(different),
        sphere_volume_relative_error=(info['geometry']['sliced_particle_volume_um3']/
                                      info['geometry']['exact_particle_volume_um3']-1))
    return dict(info=info, delta=delta, summary=summary)


def compare_cases(coarse, fine, parameter):
    same_settings(coarse['info'], fine['info'])
    a, b = coarse['info']['scene'], fine['info']['scene']
    different = [key for key in a if a[key] != b[key]]
    if parameter not in NUMERICAL or different != [parameter]:
        raise ValueError('A convergence step must change exactly one declared numerical parameter: '+str(different))
    if a[parameter] is None or b[parameter] <= a[parameter]:
        raise ValueError('Fine numerical parameter must exceed the coarse value')
    return comparison_metrics(coarse, fine, parameter, a[parameter], b[parameter])


def compare_joint_cases(coarse, fine):
    same_settings(coarse['info'], fine['info'])
    a, b = coarse['info']['scene'], fine['info']['scene']
    changed = [k for k in a if a[k] != b[k]]
    if len(changed) < 2 or any(k not in NUMERICAL for k in changed):
        raise ValueError('Joint refinement requires multiple numerical changes and no physical changes')
    if any(a[k] is None or b[k] <= a[k] for k in changed):
        raise ValueError('Every changed joint-refinement parameter must increase')
    result = comparison_metrics(coarse, fine, 'joint',
        '/'.join(str(a[k]) for k in NUMERICAL), '/'.join(str(b[k]) for k in NUMERICAL))
    result['changed_numerical_parameters'] = changed
    result['caveat'] = 'Joint step sensitivity, not separate error estimates or an error bound.'
    return result


def comparison_metrics(coarse, fine, parameter, coarse_value, fine_value):
    b = fine['info']['scene']
    result = field_change(coarse['delta'], fine['delta'], b['incidence_deg'])
    changes = []
    for old, new in zip(coarse['summary']['channels'], fine['summary']['channels']):
        x, y = old['difference_field_efficiency'], new['difference_field_efficiency']
        changes.append(x/y-1 if y else None)
    change = result.get('relative_power_weighted_complex_field_change')
    result.update(parameter=parameter, coarse_value=coarse_value, fine_value=fine_value,
        channel_difference_field_power_relative_changes=changes,
        meets_local_tolerance=(result['comparable'] and change is not None and change < .02
            and all(c is not None and abs(c) < .05 for c in changes)),
        tolerance=dict(complex_field=.02, each_channel_difference_field_power=.05))
    return result


def summarize(config):
    cases, records = {}, {}
    oracle = ROOT/config['nist_oracle'] if config.get('nist_oracle') else None
    if oracle and sha256(oracle) != config['nist_oracle_sha256']:
        raise ValueError('NIST oracle hash mismatch')
    for name, paths in config['cases'].items():
        reference_folder, sample_folder = [ROOT/paths[key] for key in ('reference', 'sample')]
        reference_info, reference = load_checked(reference_folder)
        info, sample = load_checked(sample_folder)
        cases[name] = prepare_pair(reference_info, reference, info, sample)
        usage = sample_folder.with_suffix('.resources.json')
        records[name] = dict(cases[name]['summary'], sources=paths,
            raw_npz_sha256=info['raw_npz_sha256'], input_sha256=info['input_sha256'],
            executable_sha256=info['executable_sha256'],
            resources=json.loads(usage.read_text(encoding='utf-8')) if usage.exists() else None)
        if oracle:
            records[name]['nist_background_comparison'] = nist_comparison(
                Scene(**reference_info['scene']), reference, oracle)
    comparisons = []
    for step in config['comparisons']:
        a, b = cases[step['coarse']], cases[step['fine']]
        result = compare_joint_cases(a, b) if step['parameter'] == 'joint' else compare_cases(a, b, step['parameter'])
        comparisons.append(dict(step, **result))
    return dict(status='EXPERIMENTAL_CONVERGENCE_STUDY', cases=records, comparisons=comparisons,
        nist_oracle=dict(path=config['nist_oracle'], sha256=sha256(oracle)) if oracle else None,
        limitations=[
            'Local tolerance is a step-size test, not an error bound or convergence certificate.',
            'Tests at different fixed cutoffs do not certify a joint high-cutoff, fine-slicing solution.',
            'Supercell separation and isolated-defect limits have not been checked.',
            'Plane-wave, periodic-array efficiencies are not continuous isolated-particle DSC.',
            'Material constants are controlled comparison inputs, not calibrated at 266 nm.',
            'No camera image or quantitative SNR is certified.'])


def write_markdown(output, result):
    rows = ['# S4 controlled convergence study', '', 'Status: '+result['status'], '',
        '## Declared refinement steps', '',
        '| Parameter | Coarse -> fine | Complex difference-field change | CH1 / CH2 / CH3 power changes | Local tolerance |',
        '|---|---|---:|---|---|']
    for c in result['comparisons']:
        value = c.get('relative_power_weighted_complex_field_change')
        field = '{:.5g}%'.format(100*value) if value is not None else 'Not comparable'
        powers = ' / '.join('{:.5g}%'.format(100*x) if x is not None else 'N/A'
                             for x in c['channel_difference_field_power_relative_changes'])
        rows.append('| {} | {} -> {} | {} | {} | {} |'.format(c['parameter'],
            c['coarse'], c['fine'], field, powers, 'Pass' if c['meets_local_tolerance'] else 'Fail'))
    rows += ['', 'Power changes = coarse/fine - 1. Complex fields keep their original amplitude and phase.',
        '', '## Channel Results', '',
        '| Case | CH | Sample | Reference | Difference field | Actual sample-reference |',
        '|---|---:|---:|---:|---:|---:|']
    for name, case in result['cases'].items():
        for i, ch in enumerate(case['channels'], 1):
            rows.append('| {} | {} | {:.8g} | {:.8g} | {:.8g} | {:.8g} |'.format(name, i,
                ch['sample_efficiency'], ch['reference_efficiency'],
                ch['difference_field_efficiency'], ch['actual_power_difference']))
    rows += ['', 'All efficiencies use incident supercell power; NA=0.6, tilts=(+45,-45,0) deg.',
        'Difference-field power is |E_sample-E_reference|^2, NOT sample power minus reference power.',
        '', '## Limitations', '']+['- '+x for x in result['limitations']]
    (output/'CONVERGENCE.md').write_text('\n'.join(rows)+'\n', encoding='utf-8')


def plot(output, result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    steps = [c for c in result['comparisons'] if c.get('relative_power_weighted_complex_field_change')
             is not None and all(x is not None for x in c['channel_difference_field_power_relative_changes'])]
    if not steps:
        return
    labels = ['{}: {} to {}'.format(c['parameter'], c['coarse_value'], c['fine_value']) for c in steps]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout='constrained')
    x = np.arange(len(steps))
    heights = [100*c['relative_power_weighted_complex_field_change'] for c in steps]
    axes[0].bar(x, heights, color='#007D8A', width=.55)
    axes[0].axhline(2., color='#333333', linestyle='--', label='Local test: 2%')
    for i, value in enumerate(heights):
        axes[0].annotate('{:.3f}%'.format(value), (i, value), xytext=(0, 5),
                         textcoords='offset points', ha='center')
    axes[0].set(xticks=x, xticklabels=labels, ylabel='Weighted complex-field change (%)',
                title='Defect difference field (no phase or amplitude fit)')
    axes[0].margins(y=.25)
    axes[0].legend(frameon=False)
    for i, (label, color) in enumerate(zip(('CH1', 'CH2', 'CH3'), ('#007D8A', '#CE6939', '#535353'))):
        axes[1].bar(x+(i-1)*.23, [100*c['channel_difference_field_power_relative_changes'][i]
                    for c in steps], .23, label=label, color=color)
    axes[1].axhline(0, color='#333333', linewidth=.7)
    axes[1].axhline(5, color='#777777', linestyle='--', linewidth=.7)
    axes[1].axhline(-5, color='#777777', linestyle='--', linewidth=.7)
    axes[1].set(xticks=x, xticklabels=labels, ylabel='Power change: coarse/fine - 1 (%)',
                title='Three independent NA=0.6 channels')
    axes[1].legend(frameon=False)
    for axis in axes:
        axis.tick_params(axis='x', labelsize=8)
        axis.grid(axis='y', alpha=.2)
        axis.set_axisbelow(True)
    fig.savefig(output/'convergence.png', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding='utf-8'))
    result = summarize(config)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'convergence.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    write_markdown(args.output, result)
    if args.plot:
        plot(args.output, result)
    print(json.dumps(dict(status=result['status'], comparisons=result['comparisons']), indent=2))


if __name__ == '__main__':
    main()
