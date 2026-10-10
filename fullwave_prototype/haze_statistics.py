"""Solid-angle aperture integrals and uncertainty of independent intensities."""
import numpy as np
from scipy.stats import t

from .run import hemisphere_integral


def confidence_interval(values, confidence=.95):
    values = np.asarray(values, dtype=float)
    if (values.ndim != 1 or not len(values) or not np.isfinite(values).all() or
            not np.isfinite(confidence) or not 0 < confidence < 1):
        raise ValueError('Statistics require a finite, nonempty sample vector')
    n = len(values)
    mean = float(values.mean())
    if n < 2:
        return dict(samples=n, mean_um2=mean, standard_error_um2=None,
                    ci95_um2=None, relative_ci95_half_width=None)
    standard_error = float(values.std(ddof=1)/np.sqrt(n))
    half = float(t.ppf((1+confidence)/2, n-1)*standard_error)
    return dict(samples=n, mean_um2=mean, standard_error_um2=standard_error,
                ci95_um2=[mean-half, mean+half],
                relative_ci95_half_width=(half/abs(mean) if mean != 0 else
                                          0. if half == 0 else None))


def observable_samples(maps, theta, phi, apertures=(.5, .5, .5), tilts=(45, -45, 0)):
    values = np.asarray(maps, dtype=float)
    theta, phi = np.asarray(theta), np.asarray(phi)
    if (values.ndim != 3 or values.shape[1:] != (len(theta), len(phi)) or
            len(phi) < 4 or not np.isfinite(values).all() or np.any(values < 0)):
        raise ValueError('Observables require finite nonnegative intensity maps on one angular grid')
    if not np.allclose(np.diff(phi), 2*np.pi/len(phi), atol=1e-10):
        raise ValueError('Azimuth samples must be uniformly spaced over 2 pi, without endpoint duplication')
    if (len(theta) < 2 or np.any(np.diff(theta) <= 0) or
            not np.isclose(theta[0], 0) or not np.isclose(theta[-1], np.pi/2)):
        raise ValueError('Polar samples must span the upper hemisphere from 0 to pi/2')
    if (len(apertures) != 3 or len(tilts) != 3 or
            any(not np.isfinite(v) or not 0 <= v < 1 for v in apertures) or
            any(not np.isfinite(v) or abs(v) > 89.9 for v in tilts)):
        raise ValueError('Three valid NA [0,1) and tilt [-89.9,89.9] values required')
    from .aperture_integrals import integrate_maps
    caps = integrate_maps(values, theta, phi, apertures, tilts)
    hemispheres = np.array([hemisphere_integral(value, theta) for value in values])
    return np.column_stack((hemispheres, caps))


def statistical_report(maps, theta, phi, scale=1., apertures=(.5, .5, .5),
                       tilts=(45, -45, 0), target_relative_ci=.05):
    if not np.isfinite(target_relative_ci) or not 0 < target_relative_ci < 1:
        raise ValueError('Target relative confidence half-width must be in (0,1)')
    series = observable_samples(maps, theta, phi, apertures, tilts)*scale
    return report_from_observables(series, apertures, tilts, target_relative_ci)


def report_from_observables(series, apertures=(.5, .5, .5), tilts=(45, -45, 0),
                            target_relative_ci=.05):
    """Same uncertainty calculation for independently integrated field samples."""
    series = np.asarray(series, dtype=float)
    if series.ndim != 2 or series.shape[1] != 4 or not len(series) or not np.isfinite(series).all():
        raise ValueError('Four finite observables per independent sample required')
    if not np.isfinite(target_relative_ci) or not 0 < target_relative_ci < 1:
        raise ValueError('Target confidence half-width must be in (0,1)')
    names = ('hemisphere', 'ch1', 'ch2', 'ch3')
    checkpoints = sorted(set([n for n in (2, 4, 8, 16, 32, 64, 128, 256)
                              if n <= len(series)]+[len(series)]))
    summaries = []
    for count in checkpoints:
        summaries.append(dict(samples=count, observables={name: confidence_interval(series[:count, i])
                                                          for i, name in enumerate(names)}))
    final = summaries[-1]['observables']
    active = ['hemisphere']+['ch'+str(i+1) for i, na in enumerate(apertures) if na > 0]
    converged = len(series) >= 3 and all(
        final[name]['relative_ci95_half_width'] is not None and
        final[name]['relative_ci95_half_width'] <= target_relative_ci for name in active)
    simultaneous = {name: confidence_interval(series[:, names.index(name)], 1-.05/len(active)) for name in active}
    familywise_met = len(series) >= 3 and all(item['relative_ci95_half_width'] is not None and
        item['relative_ci95_half_width'] <= target_relative_ci for item in simultaneous.values())
    return dict(observables=final, prefix_statistics=summaries,
                sample_observables_um2=series.tolist(), numerical_apertures=list(apertures),
                tilt_degrees=list(tilts), target_relative_ci95_half_width=target_relative_ci,
                statistical_target_met=bool(converged),
                simultaneous_observables=simultaneous, simultaneous_statistical_target_met=bool(familywise_met),
                simultaneous_interval_method='Bonferroni familywise 95% Student t intervals over active observables; '
                                             'finite-sample approximation for independent surface realizations',
                aperture_integration='Cone-aligned Gauss solid-angle quadrature, order 20; '
                                     'periodic bilinear DSC interpolation; upper hemisphere only',
                interval_method='Student t, pointwise 95%; finite-sample approximation; '
                                'assumes independent surfaces; no simultaneous channel guarantee')


def estimated_sample_budget(report, target=.05):
    """Pilot-CV planning only, not a promise or an optional-stopping confidence rule."""
    if not np.isfinite(target) or not 0 < target < 1:
        raise ValueError('Budget target must be in (0,1)')
    estimates = {}
    for name, item in report['observables'].items():
        n, mean, se = item['samples'], abs(item['mean_um2']), item['standard_error_um2']
        if n < 3 or se is None or mean == 0:
            estimates[name] = None
            continue
        cv = se*np.sqrt(n)/mean
        total = max(3, int(np.ceil((1.96*cv/target)**2)))
        for _ in range(10):
            refined = max(3, int(np.ceil((t.ppf(.975, total-1)*cv/target)**2)))
            if refined == total:
                break
            total = refined
        estimates[name] = total
    active = ['hemisphere']+['ch'+str(i+1) for i, na in enumerate(report['numerical_apertures']) if na > 0]
    totals = [estimates[name] for name in active]
    needed = max(totals) if all(value is not None for value in totals) else None
    return dict(target_relative_ci95_half_width=target, per_observable=estimates,
                estimated_total_surfaces=needed,
                fits_128_surface_job=needed <= 128 if needed is not None else None,
                scope='Pilot coefficient-of-variation estimate with future Student t factor; '
                      'not a convergence guarantee; no repeated-peeking stopping rule')


def plot_statistics(path, report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout='constrained')
    for ax, name in zip(axes.ravel(), ('hemisphere', 'ch1', 'ch2', 'ch3')):
        items = report['prefix_statistics']
        means = [item['observables'][name]['mean_um2'] for item in items]
        errors = [(item['observables'][name]['ci95_um2'][1]-means[i])
                  if item['observables'][name]['ci95_um2'] is not None else 0
                  for i, item in enumerate(items)]
        ax.errorbar([item['samples'] for item in items], means, yerr=errors,
                    marker='o', capsize=4, color='tab:green')
        ax.set(title=name, xlabel='Independent surfaces', ylabel='DSC integral (um2)')
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
        ax.grid(alpha=.25)
    fig.suptitle('Haze sample convergence | pointwise 95% Student t intervals')
    fig.savefig(path, dpi=150)
    plt.close(fig)
