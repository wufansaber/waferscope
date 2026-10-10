"""Absolute-scale comparison of independent FEM haze and NIST PSD."""
from pathlib import Path

import numpy as np

from .angular_plotting import angular_disk_edges, format_angular_disk
from .plot_directions import mark_illumination_source_xy


def plot_ensemble(path, theta, phi, mean, standard_error, reference, report,
                  finite_reference=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import PowerNorm

    path = Path(path)
    x_edges, y_edges = angular_disk_edges(theta, phi)
    columns = 2+int(reference is not None)+int(finite_reference is not None)
    fig, axes = plt.subplots(1, columns, figsize=(5.3*columns, 5.8),
                             layout='constrained', squeeze=False)
    axes = axes[0]
    metrics = report['metrics']
    maximum = max(float(np.max(mean)),
                  float(np.max(reference)) if reference is not None else 0.,
                  float(np.max(finite_reference)) if finite_reference is not None else 0., 1e-30)
    panels = [('fem', mean, 'FEM haze | {} surfaces\nHemisphere DSC = {:.4g} um^2'.format(
        metrics['samples'], metrics['fem_mean_hemisphere_dsc_um2']))]
    if finite_reference is not None:
        panels.append(('finite', finite_reference,
                       'Finite FEM patch | expected PSD\nHemisphere DSC = {:.4g} um^2'.format(
                           metrics['finite_patch_hemisphere_dsc_um2'])))
    if reference is not None:
        panels.append(('nist', reference, 'NIST | infinite {} PSD\nHemisphere DSC = {:.4g} um^2'.format(
            metrics.get('psd_model', 'Gaussian').upper(), metrics['nist_hemisphere_dsc_um2'])))
    panels.append(('stderr', standard_error if metrics['samples'] > 1 else np.zeros_like(mean),
                   'FEM standard error | {} surfaces'.format(metrics['samples'])))
    for ax, (name, values, title) in zip(axes, panels):
        local_max = (maximum if name != 'stderr' else
                     max(float(np.nanmax(values)), 1e-30))
        image = ax.pcolormesh(x_edges, y_edges, values,
                              cmap='viridis', norm=PowerNorm(.5, vmin=0, vmax=local_max),
                              shading='flat', rasterized=True)
        format_angular_disk(ax)
        mark_illumination_source_xy(ax, report['parameters']['incidence_deg'])
        ax.set_title(title, fontsize=11)
        fig.colorbar(image, ax=ax, shrink=.72, label='DSC (um^2 / sr)')
    fig.suptitle('Haze only | global kx/k, ky/k | equivalent area {:.4g} um^2{}'.format(
        metrics['comparison_area_um2'],
        ' | shared absolute DSC scale' if reference is not None else ''), fontsize=12)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    for name, values, title in panels:
        fig, ax = plt.subplots(figsize=(7, 6.7), layout='constrained')
        local_max = maximum if name != 'stderr' else max(float(np.nanmax(values)), 1e-30)
        image = ax.pcolormesh(x_edges, y_edges, values, cmap='viridis',
                              norm=PowerNorm(.5, vmin=0, vmax=local_max),
                              shading='flat', rasterized=True)
        format_angular_disk(ax)
        mark_illumination_source_xy(ax, report['parameters']['incidence_deg'])
        ax.set_title(title, fontsize=11)
        fig.colorbar(image, ax=ax, shrink=.75, label='DSC (um^2 / sr)')
        fig.savefig(path.with_name('haze_ensemble_'+name+'.png'), dpi=150)
        plt.close(fig)
