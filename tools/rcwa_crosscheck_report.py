"""Matched analytic-geometry S4/grcwa comparison, with no fitted amplitude or phase."""
import argparse
import json
from pathlib import Path

import numpy as np

from rcwa_prototype.comparison import compare_pair, field_change
from rcwa_prototype.geometry import Scene
from s4_prototype.solver import ROOT, sha256
from tools.s4_validation_report import load


def checked_load(folder):
    info, data = load(folder)
    if sha256(folder/'modes.npz') != info['raw_npz_sha256']:
        raise ValueError('Raw NPZ hash mismatch: '+str(folder))
    if info['backend'] == 'S4' and sha256(folder/'input.lua') != info['input_sha256']:
        raise ValueError('S4 input hash mismatch')
    if not np.all(np.isfinite(data['electric_reflected_z0'])):
        raise ValueError('Nonfinite reflected field')
    return info, data


def check_methods(s4, grcwa):
    if s4['backend'] != 'S4' or grcwa['backend'] != 'grcwa':
        raise ValueError('Expected S4 and grcwa, not a relabeled backend')
    unused = ('grid_x', 'grid_y', 'grating_method')
    if any(s4['scene'][k] != grcwa['scene'][k] for k in s4['scene'] if k not in unused):
        raise ValueError('Matched physical inputs, slicing and retained orders required')
    if s4['polarization'] != 'p' or s4['factorization'] != 'default' or s4['one_dimensional']:
        raise ValueError('Crosscheck currently supports 2D p/default S4 only')
    if s4.get('epsilon_method', 'analytic') != 'analytic' or not s4.get('stripe_rule'):
        raise ValueError('S4 analytic geometry and stripe adapter required')
    if grcwa.get('particle_method') != 'analytic' or grcwa['scene']['grating_method'] != 'analytic_li':
        raise ValueError('grcwa must use matched analytic spheres and straight-grating rule')
    if grcwa['scene']['factorization'] != 'laurent':
        raise ValueError('grcwa Laurent rule required')


def field_difference(other, reference, incidence):
    result = field_change(other, reference, incidence)
    if result['comparable']:
        a = {tuple(g): i for i, g in enumerate(other['orders']) if other['propagating'][i]}
        b = {tuple(g): i for i, g in enumerate(reference['orders']) if reference['propagating'][i]}
        keys = sorted(a)
        np.testing.assert_allclose(other['directions'][[a[k] for k in keys]],
                                   reference['directions'][[b[k] for k in keys]], atol=1e-12, rtol=1e-12)
    return dict(comparable=result['comparable'],
        common_propagating_orders=result.get('common_propagating_orders'),
        power_weighted_relative_complex_field_difference=
            result.get('relative_power_weighted_complex_field_change'),
        denominator='S4 field norm; no fitted amplitude or phase')


def summarize(paths):
    loaded = {name: checked_load(ROOT/path) for name, path in paths.items()}
    for name in ('sample', 'reference'):
        check_methods(loaded['s4_'+name][0], loaded['grcwa_'+name][0])
    data = {}
    for backend in ('s4', 'grcwa'):
        sample_info, sample = loaded[backend+'_sample']
        reference_info, reference = loaded[backend+'_reference']
        if sample_info['scene'] != reference_info['scene']:
            raise ValueError('Matched sample/reference inputs required')
        if not sample_info['particle'] or reference_info['particle']:
            raise ValueError('Sample/reference particle flags are invalid')
        scene = Scene(**sample_info['scene'])
        data[backend] = dict(summary=compare_pair(reference, sample, scene),
            sample=sample, reference=reference,
            delta=dict(sample, electric_reflected_z0=
                       sample['electric_reflected_z0']-reference['electric_reflected_z0']))
    differences = {name: field_difference(data['grcwa'][name], data['s4'][name], scene.incidence_deg)
                   for name in ('sample', 'reference', 'delta')}
    return dict(status='SAME_CUTOFF_IMPLEMENTATION_CROSSCHECK', fields=differences,
        channels={key: value['summary'] for key, value in data.items()},
        sources={key: dict(path=path, npz_sha256=loaded[key][0]['raw_npz_sha256'])
                 for key, path in paths.items()},
        limitations=['Both are Fourier modal methods; this is not an independent physical-method certification.',
            'Agreement at one cutoff does not certify slicing, cell separation or a true isolated defect.',
            'The remaining complex-field difference is reported, not fitted away.'])


def write_outputs(output, result):
    rows = ['# S4/grcwa implementation crosscheck', '',
        '| Field | Weighted relative complex-field difference |', '|---|---:|']
    for name, data in result['fields'].items():
        value = data['power_weighted_relative_complex_field_difference']
        rows.append('| {} | {} |'.format(name, '{:.7g}%'.format(value*100) if value is not None else 'N/A'))
    rows += ['', 'The denominator is the original S4 field norm. No phase or amplitude fit.', '',
        '| Backend | CH | Sample | Reference | Difference field | Actual power difference |',
        '|---|---:|---:|---:|---:|---:|']
    for backend, data in result['channels'].items():
        for i, ch in enumerate(data['channels'], 1):
            rows.append('| {} | {} | {:.9g} | {:.9g} | {:.9g} | {:.9g} |'.format(backend, i,
                ch['sample_efficiency'], ch['reference_efficiency'],
                ch['difference_field_efficiency'], ch['actual_power_difference']))
    rows += ['', 'All powers are efficiencies relative to the same incident supercell power.', '']
    rows += ['- '+text for text in result['limitations']]
    (output/'CROSSCHECK.md').write_text('\n'.join(rows)+'\n', encoding='utf-8')
    plot(output, result)


def plot(output, result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout='constrained')
    names = [name for name, row in result['fields'].items()
             if row['power_weighted_relative_complex_field_difference'] is not None]
    values = [100*result['fields'][name]['power_weighted_relative_complex_field_difference'] for name in names]
    axes[0].bar(names, values, color='#007D8A', width=.55)
    for i, value in enumerate(values):
        axes[0].annotate('{:.3g}%'.format(value), (i, value), xytext=(0, 5),
                         textcoords='offset points', ha='center')
    axes[0].set(ylabel='Complex-field difference (%)', title='Same retained modes and analytic geometry')
    axes[0].margins(y=.3)
    x = np.arange(3)
    for shift, key, color in ((-.17, 's4', '#007D8A'), (.17, 'grcwa', '#CE6939')):
        axes[1].bar(x+shift, [c['difference_field_efficiency'] for c in result['channels'][key]['channels']],
                     .34, label=key, color=color)
    axes[1].set(xticks=x, xticklabels=['CH1', 'CH2', 'CH3'], ylabel='Difference-field collection efficiency',
                title='Independent NA=0.6 collection channels')
    axes[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    axes[1].legend(frameon=False)
    for ax in axes:
        ax.grid(axis='y', alpha=.2)
        ax.set_axisbelow(True)
    fig.savefig(output/'crosscheck.png', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    paths = json.loads(args.config.read_text(encoding='utf-8'))
    result = summarize(paths)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'crosscheck.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    write_outputs(args.output, result)
    print(json.dumps(result['fields'], indent=2))


if __name__ == '__main__':
    main()
