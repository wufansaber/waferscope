"""CH3 raw simulated camera triptych from an explicitly audited complex source."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np


def raw_triplet(arrays):
    sample = np.asarray(arrays['detector_gray_first_sample_dn'])
    reference = np.asarray(arrays['detector_gray_first_reference_dn'])
    if sample.shape != reference.shape or sample.ndim != 2:
        raise ValueError('Matching two-dimensional raw frames required')
    return sample, reference, sample.astype(np.int32)-reference.astype(np.int32)


def render(output, arrays, report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sample, reference, difference = raw_triplet(arrays)
    detector = report['settings']['detector']
    dx = detector['pixel_pitch_um']
    nx, ny = sample.shape[1], sample.shape[0]
    extent = (-nx*dx/2, nx*dx/2, -ny*dx/2, ny*dx/2)
    lo, hi = min(sample.min(), reference.min()), max(sample.max(), reference.max())
    limit = max(float(abs(difference).max()), 1)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.6), constrained_layout=True)
    for ax, frame, name in zip(axes, (sample, reference, difference),
                                ('Sample', 'Reference', 'Sample - Reference')):
        signed = name.startswith('Sample -')
        im = ax.imshow(frame, origin='lower', extent=extent, interpolation='nearest',
            cmap='RdBu_r' if signed else 'gray', vmin=-limit if signed else lo,
            vmax=limit if signed else max(hi, lo+1))
        ax.set(title='CH3 | '+name, xlabel='Sensor x (um)', ylabel='Sensor y (um)')
        fig.colorbar(im, ax=ax, shrink=.72, label='Signed DN' if signed else 'Raw DN')
    fig.suptitle('SIMULATED camera | 100 nm particle | 40x | NA 0.6 | 266 nm\n'
                 '1 us | 5 um pixels | one raw pair, no averaging\n'
                 'FEM convergence NOT passed: diagnostic preview only', fontsize=12)
    fig.savefig(output/'CH3_sample_reference_difference.png', dpi=180)
    plt.close(fig)


def check_image_propagation(arrays, p, settings):
    from fullwave_prototype.direct_imaging import direct_image
    from tools.patterned_numerical_audit import relative
    axis = arrays['image_axis_um']
    indices = np.flatnonzero(abs(axis) <= 20)[::5]
    x, y = np.meshgrid(axis[indices], axis[indices])
    u = arrays['pupil_u_axis']
    uu, vv = np.meshgrid(u, u)
    mask = arrays['pupil_mask']
    uv = np.column_stack((uu[mask], vv[mask]))
    weights = np.full(len(uv), (u[1]-u[0])**2)
    errors = {}
    for state in ('with', 'without'):
        _, _, flux = direct_image(arrays['amplitude_'+state][mask], uv, weights,
            np.column_stack((x.ravel(), y.ravel())), p.wavelength_nm/1000,
            settings.tilt_degrees[2], settings.magnification, settings.image_defocus_um,
            settings.transmission, settings.analyzer)
        errors[state] = relative(flux.reshape(x.shape)*settings.irradiance_w_um2,
                                arrays['image_'+state+'_w_um2'][np.ix_(indices, indices)])
    return dict(flux_relative_l2=errors, passed=max(errors.values()) < 1e-10,
        scope='Independent direct-vector propagation versus FFT on the same sampled pupil; '
              'not a check of FEM, source or pupil convergence.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--pupil-samples', type=int, choices=(32, 64), default=64)
    args = parser.parse_args()
    from fullwave_prototype.parameters import Parameters
    from fullwave_prototype.imaging import ImagingSettings, calculate_channel
    from detector_model import DetectorParameters
    from fullwave_prototype.radiation import radiation_amplitude
    from defect_detection import calculate_detection
    from gray_acquisition import acquire_gray
    audit = json.loads((args.source/'audit.json').read_text(encoding='utf-8'))
    p = Parameters(**audit['parameters'])
    config = dict(audit['settings'])
    config['detector'] = DetectorParameters(**config['detector'])
    settings = replace(ImagingSettings(**config), pupil_samples=args.pupil_samples)
    settings.validate()
    config['pupil_samples'] = args.pupil_samples
    config['detector'] = audit['settings']['detector']
    sources = {}
    for name in ('with', 'without'):
        with np.load(args.source/('sources_'+name+'.npz')) as data:
            sources[name] = data['points'], data['weights']
    def provider(directions, state):
        print('CH3 {}: propagating {} volume sources to {} pupil rays'.format(
            state, len(sources[state][0]), len(directions)), flush=True)
        return radiation_amplitude(*sources[state], directions, p)
    arrays, metrics = calculate_channel({}, p, settings, 2, amplitude_provider=provider)
    detector = calculate_detection(dict(image_axis_um=arrays['image_axis_um'],
        irradiance_w_um2=arrays['image_with_w_um2'], power_w=metrics['with']['image_power_w']),
        dict(image_axis_um=arrays['image_axis_um'], irradiance_w_um2=arrays['image_without_w_um2'],
             power_w=metrics['without']['image_power_w']), p.wavelength_nm, settings.detector,
        rng=np.random.default_rng(20261019), sensor_index=2)
    detector.update(acquire_gray(detector, settings.detector, settings.gray_snr_pairs,
                                settings.gray_snr_seed, 2))
    arrays.update({'detector_'+k: v for k, v in detector.items() if isinstance(v, np.ndarray)})
    report = dict(settings=config, parameters=audit['parameters'], pattern=audit['pattern'],
        channel=3, metrics=metrics, detector={k: v for k, v in detector.items() if not isinstance(v, np.ndarray)},
        validation_status='not_converged', source=str(args.source.resolve()),
        limitations='Simulated raw camera DN, not a real acquisition. Numerical and material accuracy not certified. '
                    'Pupil sampling convergence not established; 1e-6 optical transmission is an assumed attenuation.')
    report['image_propagation_check'] = check_image_propagation(arrays, p, settings)
    args.output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output/'ch3_fields.npz', **arrays)
    (args.output/'ch3_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    render(args.output, arrays, report)
    print('CH3 diagnostic preview saved: '+str(args.output.resolve()), flush=True)


if __name__ == '__main__':
    main()
