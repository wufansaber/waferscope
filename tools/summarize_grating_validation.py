"""Evidence report for periodic benchmarks and finite-surface consistency."""
import argparse
import json
from pathlib import Path

import numpy as np


def relative(a, b):
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Non-finite comparison')
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b), 1e-300))


def load_arrays(path):
    with np.load(path, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}


def metrics(actual, reference, efficiency='efficiency', amplitude='amplitude'):
    mask = reference['direction'][:, 2] > 0
    np.testing.assert_array_equal(actual['orders'], reference['orders'])
    ray = reference['direction'][mask]
    np.testing.assert_allclose(actual['tangent'][mask], ray[:, 0], atol=1e-12)
    amp = actual[amplitude]
    if amp.ndim == 2:
        amp = amp[0]
    x, y = actual[efficiency][mask], reference['efficiency'][mask]
    cap = abs(ray[:, 0]) <= .6
    return dict(orders=reference['orders'][mask].tolist(), fem_efficiency=x.tolist(),
        nist_efficiency=y.tolist(), relative_power_error=((x-y)/y).tolist(),
        absolute_efficiency_error=(x-y).tolist(),
        total_reflectance_relative_error=float(x.sum()/y.sum()-1),
        ch3_efficiency=float(x[cap].sum()), nist_ch3_efficiency=float(y[cap].sum()),
        ch3_relative_power_error=float(x[cap].sum()/y[cap].sum()-1),
        complex_relative_l2=relative(amp[mask], reference['amplitude'][mask]))


def plot(output, report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cases = [('NIST RCW', None, '#333333'), ('2D FEM', 'scalar', '#287a77'),
             ('3D FEM: plane', 'vector_plane', '#b65048'), ('3D FEM: volume', 'vector_volume', '#8a66a7')]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    base = report['scalar']
    x = np.arange(len(base['orders']))
    for i, (label, key, color) in enumerate(cases):
        value = base['nist_efficiency'] if key is None else report[key]['fem_efficiency']
        axes[0].bar(x+(i-1.5)*.19, value, width=.18, label=label, color=color)
    axes[0].set(xticks=x, xticklabels=base['orders'], xlabel='Canonical diffraction order m',
        ylabel='Reflected / incident normal power per period', yscale='log',
        ylim=(1e-4, 1), title='Periodic patterned wafer, no particle')
    axes[0].legend(fontsize=8)
    for label, key, color in cases[1:]:
        axes[1].plot(base['orders'], np.array(report[key]['relative_power_error'])*100,
                     'o-', color=color, label=label)
    axes[1].axhline(0, color='#555555', linewidth=.7)
    for value in (-5, 5):
        axes[1].axhline(value, color='#777777', linewidth=.7, linestyle='--')
    axes[1].set(xlabel='Canonical diffraction order m', ylabel='Relative order-power difference (%)',
                title='No fitted phase, amplitude scaling or image alignment')
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('266 nm | pitch 500 nm | ridge 250 nm | depth 100 nm | p, 75 deg\n'
        'Infinite periodic benchmark only; not certification of finite SEM patterns', fontsize=12)
    fig.savefig(output/'periodic_model_comparison.png', dpi=180)
    plt.close(fig)
    finite = report.get('finite_surface_consistency')
    if finite:
        fig, ax = plt.subplots(figsize=(7, 4.6), constrained_layout=True)
        for state, name in (('without', 'Reference'), ('with', 'Sample')):
            checks = finite[state]
            ax.plot([c['order'] for c in checks],
                [100*c['volume_vs_surface_relative_l2'] for c in checks], 'o-', label=name)
        ax.axhline(5, color='#777777', linestyle='--', label='5% selected-check target')
        ax.set(xlabel='Closed-surface Gauss order per interval', ylabel='Complex far-field relative L2 (%)',
            title='Finite patterned wafer: volume vs reciprocity surface\n'
                  'Same FEM field; this does not establish mesh convergence')
        maximum = max(100*c['volume_vs_surface_relative_l2']
                      for rows in finite.values() for c in rows)
        ax.set_ylim(0, max(10, maximum*1.15))
        ax.legend()
        fig.savefig(output/'finite_near_to_far_consistency.png', dpi=180)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--vector', default='vector_p3_densepml')
    args = parser.parse_args()
    root = args.root
    reference = load_arrays(root/'nist245.npz')
    old_ref = load_arrays(root/'nist145.npz')
    scalar = load_arrays(root/'explicit_p4_refine'/'modes.npz')
    vector = load_arrays(root/args.vector/'modes.npz')
    scalar_ref = dict(scalar, direction=reference['direction'], amplitude=scalar['amplitude'][0])
    scalar_mesh = load_arrays(root/'explicit_p4'/'modes.npz')
    scalar_domain = load_arrays(root/'explicit_p4_domain'/'modes.npz')
    report = dict(scope='Numerical verification with unchanged illustrative optical constants, not experimental calibration.',
        inputs=json.loads((root/'nist245.json').read_text(encoding='utf-8'))['inputs'],
        reference_convergence=metrics(dict(old_ref, tangent=scalar['tangent']), reference),
        scalar=metrics(scalar, reference), vector_plane=metrics(vector, reference),
        vector_volume=metrics(vector, reference, 'volume_efficiency_q11', 'volume_amplitude_q11'),
        scalar_checks=json.loads((root/'explicit_p4_refine'/'report.json').read_text(encoding='utf-8'))['physics'],
        vector_checks=json.loads((root/args.vector/'report.json').read_text(encoding='utf-8'))['checks'],
        provenance=dict(root=str(root.resolve()), scalar='explicit_p4_refine', vector=args.vector,
                        nist_fourier_orders=[145, 245]),
        scalar_mesh_sensitivity=metrics(scalar_mesh, scalar_ref),
        scalar_domain_and_remesh_sensitivity=metrics(scalar_domain, scalar_ref),
        normalization='Efficiencies per infinite periodic cell, not finite-particle DSC in um^2. '
            'CH3 is NA 0.6 about the normal and includes only propagating orders within this cone.',
        coordinate_mapping='NIST grating transverse +x maps to application global +KY; canonical m=-NIST index.',
        independent_references=['https://pages.nist.gov/SCATMECH/docs/rcwmodel.htm',
            'https://pages.nist.gov/pySCATMECH/UsingRCW.html'],
        limitations=['Infinite periodic lines without a particle are not the previous finite two-trench scene.',
            'No reconstruction or radiometric fitting to the supplied SEM/dark-field screenshots.',
            'Legacy optical constants remain 5.57+0.387i at 266 nm solely for controlled numerical comparison.',
            'Camera images, finite defects and phase-sensitive particle/pattern coupling remain unverified.'])
    report['status'] = dict(scalar_periodic='passed_selected_checks' if
        max(abs(v) for v in report['scalar']['relative_power_error']) < .005 and
        abs(report['scalar_checks']['energy_closure_error']) < .001 else 'failed',
        vector_periodic='passed_selected_5pct_order_power_checks' if
        max(abs(v) for v in report['vector_volume']['relative_power_error']) < .05 and
        max(abs(v) for v in report['vector_plane']['relative_power_error']) < .05 else 'not_all_orders_within_5pct',
        finite_pattern='not_certified')
    finite_path = root/'finite_surface_p2'/'audit.json'
    if finite_path.exists():
        finite = json.loads(finite_path.read_text(encoding='utf-8'))
        report['finite_surface_consistency'] = {k: v['reciprocity_surface_checks']
            for k, v in finite['rows'][-1]['states'].items()}
        report['finite_parameters'] = finite['parameters']
        report['finite_geometry'] = finite['pattern']
        report['finite_interpretation'] = ('This is a same-field consistency discrepancy, not an absolute FEM error bound. '
            'Surface quadrature is not fully converged; inspect its 48-to-72 change separately. '
            'Mesh, corner refinement, enclosure position and PML convergence remain required.')
    resource = root/'vector_p3_corner.resources.json'
    if resource.exists():
        report['corner_refinement_resource_attempt'] = json.loads(resource.read_text(encoding='utf-8'))
        report['corner_refinement_note'] = ('Stopped at the configured 8 GiB private-commit accounting limit, '
            'not an observed operating-system out-of-memory event. No result from this run is used as validation evidence.')
    report['next_checks'] = ['Converge the closed-surface quadrature and repeat with another enclosure.',
        'Refine physical substrate and trench corners with a solver/memory strategy that completes safely.',
        'Require finite sample/reference mesh and PML convergence before quantitative camera or SNR claims.',
        'Calibrate 266 nm optical constants and SEM dimensions before experimental comparison.']
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'validation_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    plot(args.output, report)
    print(json.dumps(dict(status=report['status'], scalar=report['scalar'],
        vector=report['vector_volume'], output=str(args.output.resolve())), indent=2))


if __name__ == '__main__':
    main()
