"""Read measured theory audits and plot results without hiding failed cases."""
import argparse
import json
from pathlib import Path

import numpy as np


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def collect(root):
    names = ('o2_b0', 'o2_b4', 'o2_a1', 'o3_a1', 'o3_compact', 'o3_compact_loss', 'o3_uv_small')
    report = dict(sphere=[], incomplete=[], optics=read(root/'optics_fixed/report.json'),
                  bv_convergence=read(root/'bv_convergence/report.json'), flat=[])
    for name in names:
        path = root/name/'report.json'
        if path.exists():
            case = dict(case=name, **read(path))
            if (root/name/'huygens.json').exists():
                case['huygens'] = read(root/name/'huygens.json')
            report['sphere'].append(case)
        else:
            report['incomplete'].append(name)
    for pol in ('p', 's'):
        path = root/'flat_compact'/('flat_'+pol)/'report.json'
        if path.exists():
            report['flat'].append(dict(polarization=pol, **read(path)))
        else:
            report['incomplete'].append('flat_'+pol)
    path = root/'rough_balanced/imaging_report.json'
    if path.exists():
        report['rough_imaging'] = read(path)
        report['rough_resources'] = read(root/'resources_rough_balanced.json')
    else:
        report['incomplete'].append('rough_balanced')
    path = root/'rough_balanced128/imaging_report.json'
    if path.exists() and 'rough_imaging' in report:
        fine = read(path)
        coarse = report['rough_imaging']
        if fine['scattering_parameters'] != coarse['scattering_parameters'] or fine['roughness'] != coarse['roughness']:
            raise ValueError('Pupil convergence requires identical scattering parameters')
        for key, value in coarse['settings'].items():
            if key != 'pupil_samples' and fine['settings'][key] != value:
                raise ValueError('Pupil convergence settings differ: '+key)
        report['rough_imaging128'] = fine
        report['rough_pupil_convergence'] = [dict(
            tilt_deg=b['tilt_deg'], power64_pw=a['with']['image_power_w']*1e12,
            power128_pw=b['with']['image_power_w']*1e12,
            relative_power_change=abs(a['with']['image_power_w']/b['with']['image_power_w']-1))
            for a, b in zip(coarse['channels'], fine['channels'])]
        report['rough128_resources'] = read(root/'resources_rough128.json')
    else:
        report['incomplete'].append('rough_balanced128')
    path = root/'rough_convergence.json'
    if path.exists():
        report['rough_discretization'] = read(path)
    else:
        report['incomplete'].append('rough_convergence')
    report['scope'] = 'Finite benchmark set; no universal accuracy or experimental-validation claim.'
    (root/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    return report


def plot(report, root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    cases = {r['case']: r for r in report['sphere']}
    colors = ('#69757c', '#197c78')
    for i, (name, label) in enumerate((('o2_b0', 'Previous order 2'), ('o3_compact', 'Balanced order 3'))):
        data = cases[name]
        axes[0, 0].bar(np.arange(3)+(i-.5)*.35,
            [100*data[k] for k in ('scattering_relative_error', 'extinction_relative_error', 'energy_balance_relative')],
            width=.35, color=colors[i], label=label)
        axes[0, 1].plot([n['radius_nm'] for n in data['near_shells']],
            [100*n['electric_relative_l2'] for n in data['near_shells']], 'o-', color=colors[i], label=label+' E')
        axes[0, 1].plot([n['radius_nm'] for n in data['near_shells']],
            [100*n['magnetic_relative_l2'] for n in data['near_shells']], 's--', color=colors[i], label=label+' H')
    axes[0, 0].set(title='Free sphere: exact Mie and energy balance', xticks=range(3),
                   xticklabels=['Scattering', 'Extinction', 'Energy closure'], ylabel='Relative discrepancy (%)')
    axes[0, 0].legend()
    axes[0, 1].set(title='Independent complex near fields (no fitting)', xlabel='Distance from sphere center (nm)',
                   ylabel='Discrete relative L2 (%)')
    axes[0, 1].legend(fontsize=8)
    with np.load(root/'optics_fixed/optics.npz') as data:
        key = 'm1.0_z0.0'
        axes[1, 0].plot(data[key+'_axis'], data[key+'_reference'], color='#af5639', label='Independent Bessel integral')
        axes[1, 0].plot(data[key+'_axis'], data[key+'_flux'], '.', color='#197c78', label='Corrected vector FFT')
    axes[1, 0].set(title='Image NA 0.8: absolute axial flux', xlabel="Image x (um)", ylabel='Flux / incident irradiance')
    axes[1, 0].legend(fontsize=8)
    for index, pol in enumerate(report['flat']):
        axes[1, 1].bar(np.arange(3)+(index-.5)*.35,
            [100*c['image_flux_relative_l2'] for c in pol['channels']], width=.35,
            color=('#3764a0', '#bc723b')[index], label=pol['polarization']+' incidence')
    axes[1, 1].set(title='Balanced FEM vs BV: wafer images', xticks=range(3),
                   xticklabels=['+45 deg', '-45 deg', '0 deg'], ylabel='Image relative L2 (%)')
    axes[1, 1].legend()
    for ax in axes.flat:
        ax.grid(alpha=.2)
        ax.set_axisbelow(True)
    fig.suptitle('Independent theory checks and measured numerical improvements', fontsize=16)
    fig.savefig(root/'theory_validation.png', dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', default='build/theory_audit')
    args = parser.parse_args()
    root = Path(args.root)
    report = collect(root)
    plot(report, root)
    print(json.dumps(dict(incomplete=report['incomplete'], output=str(root.resolve())), indent=2))


if __name__ == '__main__':
    main()
