"""Summarize retained numerical evidence without changing any simulation source."""
import argparse
import json
from pathlib import Path
import numpy as np
from fullwave_prototype.haze_statistics import report_from_observables


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True)
    args = parser.parse_args()
    root = Path(args.root)
    def read(relative):
        return json.loads((root/relative).read_text(encoding='utf-8'))
    ensemble = read('ensemble128/whole_haze_report.json')
    values = np.asarray(ensemble['statistics']['sample_observables_um2'])
    count = ensemble['extension']['previous_samples']
    settings = ensemble['imaging_settings']
    stats = {name: report_from_observables(series, settings['numerical_apertures'], settings['tilt_degrees'])
             for name, series in [('original', values[:count]), ('all', values), ('additional_only', values[count:])]}
    physics = read('power_check/scale_0p1_order_1/report.json')
    mesh = read('surface_mesh/surface_mesh_report.json')
    beam = read('beam_inputs/beam_definitions_report.json')
    hybrid = read('large_hybrid/hybrid_report.json')
    range_path = 'range_refined/range_report.json' if (root/'range_refined/range_report.json').is_file() else 'large_range/range_report.json'
    range_report = read(range_path)
    result = dict(statistics=stats, physics=physics, spatial_taper=mesh, beam=beam,
        range=range_report, hybrid=dict(local_change_relative_l2_error=hybrid['local_change_relative_l2_error'],
        local_fem_vs_spm_complex_relative_l2=hybrid['local_fem_vs_spm_complex_relative_l2'],
        production_whole_surface_plan=hybrid['whole_surface_plan']),
        finite_beam_workflows={name: read('haze_'+name+'/whole_haze_report.json')
            for name in ('wafer_rectangle', 'measured_wafer') if
            (root/('haze_'+name+'/whole_haze_report.json')).is_file()},
        resource_reports={p.name: json.loads(p.read_text(encoding='utf-8')) for p in root.glob('*_guard.json')},
        scope='Independent campaigns, not one identical particle realization. Statistics are haze-only SPM; '
              'low-RMS absolute-power/closure is a small haze-only FEM test. '
              'Range refinement does not promote the production source window. Higher-order global coupling remains unverified.')
    (root/'campaign_summary.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), layout='constrained')
    names = ['ch1', 'ch2', 'ch3']
    x = np.arange(3)
    for label, source, offset, key in [('Original 16', 'original', -.27, 'observables'),
        ('All 128', 'all', 0, 'observables'), ('128 simultaneous', 'all', .27, 'simultaneous_observables')]:
        axes[0].bar(x+offset, [100*stats[source][key][n]['relative_ci95_half_width'] for n in names],
                    width=.25, label=label)
    axes[0].axhline(5, color='black', linestyle='--', linewidth=1)
    axes[0].set(xticks=x, xticklabels=names, ylabel='Relative interval half-width (%)',
                title='Haze-only statistics | fixed final count')
    axes[0].legend(fontsize=8)
    powers = physics['absolute_power_check']
    for label, key, offset in [('FEM', 'hemisphere_and_channels_fem_w', -.17),
                               ('SPM', 'hemisphere_and_channels_spm_w', .17)]:
        axes[1].bar(np.arange(4)+offset, np.asarray(powers[key])*1e9, width=.32, label=label)
    axes[1].set(xticks=np.arange(4), xticklabels=['Hemi', 'ch1', 'ch2', 'ch3'],
        ylabel='Scattered power (nW)', title='0.2 nm RMS | 1 W | small haze patch')
    axes[1].legend(fontsize=8)
    cases = mesh['cases'][1:]
    axes[2].bar(np.arange(4), [100*max(c['complex_relative_l2']) for c in cases], color='tab:green')
    axes[2].set(xticks=np.arange(4), xticklabels=['Coarse', 'Half taper', 'Double taper', 'Fine'],
                yscale='log', ylabel='Maximum cap complex-field change (%)',
                title='Same-band parent | not bandwidth proof')
    axes[2].axhline(1, color='black', linestyle='--', linewidth=1)
    fig.savefig(root/'campaign_summary.png', dpi=150)
    plt.close(fig)
    print(json.dumps({key: {name: 100*item['relative_ci95_half_width'] for name, item in
                     value['simultaneous_observables'].items()} for key, value in stats.items()}, indent=2))


if __name__ == '__main__':
    main()
