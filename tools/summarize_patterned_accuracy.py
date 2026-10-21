"""Absolute-field convergence evidence and an independent flat-wafer oracle."""
import argparse
import json
from pathlib import Path

import numpy as np

from tools.patterned_numerical_audit import relative


def power(amplitude, directions, weights):
    return float(np.sum(weights/directions[:, 2]*np.sum(abs(amplitude)**2, axis=1)))


def compare(a, b, directions, weights):
    results = {}
    for state in ('with', 'without', 'change'):
        x, y = ((a['with']-a['without'], b['with']-b['without']) if state == 'change'
                else (a[state], b[state]))
        ix, iy = np.sum(abs(x)**2, axis=1), np.sum(abs(y)**2, axis=1)
        results[state] = dict(complex_relative_l2=relative(x, y), intensity_relative_l2=relative(ix, iy),
            cap_dsc_coarse_um2=power(x, directions, weights),
            cap_dsc_fine_um2=power(y, directions, weights),
            cap_power_relative_error=power(x, directions, weights)/power(y, directions, weights)-1)
    return results


def load_case(folder, bonus, order):
    with np.load(folder/'amplitudes.npz') as data:
        return ({s: data['b{}_{}_q{}'.format(bonus, s, order)] for s in ('with', 'without')},
                data['directions'], data['angular_weights'])


def check_same_physics(reference, candidate):
    keys = ('wavelength_nm', 'diameter_nm', 'gap_nm', 'particle_n', 'particle_k',
            'substrate_n', 'substrate_k', 'incidence_deg', 'polarization', 'beam_profile')
    if any(reference['parameters'][k] != candidate['parameters'][k] for k in keys):
        raise ValueError('Cannot claim convergence across different physical inputs')
    keys = ('pitch_um', 'trench_width_um', 'trench_depth_um', 'line_length_um')
    if any(reference['pattern'][k] != candidate['pattern'][k] for k in keys):
        raise ValueError('Cannot claim convergence across different trench geometry')
    if reference['settings'] != candidate['settings']:
        raise ValueError('Optical/camera inputs must remain unchanged in this audit')


def flat_reference(folder):
    from fullwave_prototype.imaging_audit import bv_reference
    from fullwave_prototype.direct_imaging import polar_cap
    report = json.loads((folder/'audit.json').read_text(encoding='utf-8'))
    with np.load(folder/'amplitudes.npz') as data:
        directions, weights = data['directions'], data['angular_weights']
        key = max((k for k in data.files if k.startswith('b') and '_with_q' in k),
                  key=lambda k: (int(k.split('_')[0][1:]), int(k.rsplit('q', 1)[1])))
        actual = data[key]
    refs = [bv_reference(report['parameters'], weights, directions, folder/('bv_'+str(order)),
             r'D:\anaconda\envs\NIST2\python.exe', lmax=order) for order in (14, 18)]
    _, fine_weights, fine_directions = polar_cap(.6, 0, 16, 64)
    fine_reference = bv_reference(report['parameters'], fine_weights, fine_directions,
        folder/'bv_cap_fine', r'D:\anaconda\envs\NIST2\python.exe', lmax=18)
    i, j = [np.sum(abs(v)**2, axis=1) for v in (actual, refs[-1])]
    return dict(intensity_relative_l2=relative(i, j),
        cap_power_relative_error=power(actual, directions, weights)/power(refs[-1], directions, weights)-1,
        bv_multipole_relative_l2=relative(refs[0], refs[-1]),
        bv_cap_quadrature_relative_error=power(refs[-1], directions, weights)/
            power(fine_reference, fine_directions, fine_weights)-1,
        scope='Flat silicon wafer and sphere only. No fitted amplitude; intensity comparison is independent of Jones global phase conventions.')


def plot_evidence(output, result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    labels = ['Flat wafer P1', 'Flat wafer P2']
    values = [100*abs(result['flat_oracle'][k]['cap_power_relative_error'])
              for k in ('flat_p1', 'flat_p2')]
    axes[0].bar(labels, values, color=('#b04949', '#338369'))
    axes[0].set(title='Independent NIST sphere benchmark', ylabel='CH3 collected-power error (%)')
    for x, value in enumerate(values):
        axes[0].text(x, value+.15, '{:.2f}%'.format(value), ha='center')
    keys = ['integration_6_vs_10', 'mesh_factor_1_vs_08', 'p1_vs_p2', 'pml_02_vs_03']
    labels = ['Integration\n6 vs 10', 'Mesh\nh vs 0.8h', 'Order\nP1 vs P2', 'PML + remesh\n0.2 vs 0.3 um']
    if 'p2_h1_vs_p3_h18' in result['comparisons']:
        keys.append('p2_h1_vs_p3_h18')
        labels.append('Order + mesh\nP2/h vs P3/1.8h')
    values = [100*max(result['comparisons'][k][s]['complex_relative_l2']
                      for s in ('with', 'without')) for k in keys]
    axes[1].bar(labels, values, color=['#338369']+['#b04949']*(len(values)-1))
    for x, value in enumerate(values):
        axes[1].text(x, value+.7, '{:.2f}%'.format(value), ha='center')
    axes[1].set(title='Patterned scene: absolute complex-field sensitivity',
                ylabel='CH3 selected-direction relative L2 (%)', ylim=(0, max(values)*1.2))
    for ax in axes:
        ax.axhline(5, color='#555555', linestyle='--', linewidth=1, label='5% selected-check target')
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(fontsize=8)
    fig.suptitle('100 nm sphere | 266 nm | NA 0.6 | No haze\n'
                 'Flat benchmark improved; patterned scene is NOT converged', fontsize=13)
    fig.savefig(output/'accuracy_checks.png', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    cases = {}
    for name in ('p1_integration', 'p1_fine', 'p2', 'p1_pml03', 'p3'):
        path = root/name/'audit.json'
        if path.exists():
            cases[name] = json.loads(path.read_text(encoding='utf-8'))
            check_same_physics(cases['p1_integration'], cases[name])
    low, directions, weights = load_case(root/'p1_integration', 0, 9)
    high, _, _ = load_case(root/'p1_integration', 10, 9)
    mid, _, _ = load_case(root/'p1_integration', 6, 9)
    p2, _, _ = load_case(root/'p2', 8, 9)
    fine, _, _ = load_case(root/'p1_fine', 10, 13)
    comparisons = {name: compare(a, b, directions, weights) for name, a, b in (
        ('underintegration_0_vs_10', low, high), ('integration_6_vs_10', mid, high),
        ('p1_vs_p2', high, p2), ('mesh_factor_1_vs_08', high, fine))}
    if (root/'p1_pml03'/'audit.json').exists():
        changed, _, _ = load_case(root/'p1_pml03', 10, 11)
        comparisons['pml_02_vs_03'] = compare(high, changed, directions, weights)
    if (root/'p3'/'audit.json').exists():
        p3, _, _ = load_case(root/'p3', 10, 13)
        comparisons['p2_h1_vs_p3_h18'] = compare(p2, p3, directions, weights)
    result = dict(status='not_converged_do_not_use_for_quantitative_camera_prediction',
        comparisons=comparisons,
        input_checks='Physical defect, trench geometry, illumination and imaging inputs unchanged across patterned probes.',
        case_parameters={k: dict(parameters=v['parameters'], pattern=v['pattern']) for k, v in cases.items()},
        flat_oracle={name: flat_reference(root/name) for name in ('flat_p1', 'flat_p2')},
        thresholds=dict(complex_relative_l2=.05, cap_power_relative_error=.05,
                        source_quadrature_relative_l2=.01),
        limitations=['48 CH3 directions are probes, not a hemisphere integral convergence test.',
            'P2 physical region uses P1 PML: its order must be checked independently.',
            'Changing PML thickness regenerates the mesh; the sensitivity is not an isolated PML reflection estimate.',
            'No fitted field phase, amplitude scaling or forced image symmetry in the comparisons.',
            'Silicon 5.57+0.387i and sphere 1.59+0i are legacy illustrative inputs, not calibrated at 266 nm.',
            'Two finite trenches do not represent an infinite grating or a reconstruction of SEM.',
            'No domain-size or energy-balance certification. Low algebraic residual is not a physical error estimate.'])
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'accuracy_report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    plot_evidence(args.output, result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
