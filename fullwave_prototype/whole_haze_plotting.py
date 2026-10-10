"""Absolute-unit whole-beam roughness, far field and mean images."""
import numpy as np


def plot_results(output, arrays, report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import PowerNorm
    from .angular_plotting import angular_disk_edges, format_angular_disk
    from .haze_statistics import plot_statistics

    axis = arrays['x_um']
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout='constrained')
    for ax, key, title, unit, cmap in (
        (axes[0], 'first_height_nm', 'First whole-beam rough surface', 'Height (nm)', 'coolwarm'),
        (axes[1], 'incident_flux_w_um2', 'Incident beam on modeled rough wafer', 'W / um2', 'viridis')):
        image = ax.pcolormesh(axis, axis, arrays[key], shading='auto', cmap=cmap)
        ax.set(xlabel='Wafer X (um)', ylabel='Wafer Y (um)', aspect='equal', title=title)
        fig.colorbar(image, ax=ax, label=unit)
    fig.suptitle('First-order surface model | taper-weighted coverage {:.3%}'.format(
        report['plan']['checked_taper_weighted_power_coverage']))
    fig.savefig(output/'whole_surface.png', dpi=140)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 6), layout='constrained')
    x, y = angular_disk_edges(arrays['theta'], arrays['phi'])
    values = arrays['mean_dsc_um2_sr']
    image = ax.pcolormesh(x, y, values, shading='flat', cmap='viridis',
                         norm=PowerNorm(.5, vmin=0, vmax=max(float(values.max()), 1e-300)))
    format_angular_disk(ax)
    ax.set_title('Whole-beam haze | first-order | {} surfaces\nHemisphere mean power {:.4g} W'.format(
        report['metrics']['samples'], report['metrics']['mean_power_w']['hemisphere']))
    ax.scatter([0], [-np.sin(np.deg2rad(report['parameters']['incidence_deg']))],
               marker='*', c='white', edgecolor='black', zorder=8)
    fig.colorbar(image, ax=ax, label='DSC (um2 / sr)')
    fig.savefig(output/'whole_haze_far.png', dpi=140)
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), layout='constrained')
    vmax = max(float(arrays['ch{}_mean_image_w_um2'.format(i)].max()) for i in (1, 2, 3))
    for i, ax in enumerate(axes, 1):
        image_axis = arrays['ch{}_image_axis_um'.format(i)]
        image = ax.pcolormesh(image_axis, image_axis, arrays['ch{}_mean_image_w_um2'.format(i)],
            cmap='inferno', shading='auto', norm=PowerNorm(.5, vmin=0, vmax=max(vmax, 1e-300)))
        ax.set(xlabel='Image X (um)', ylabel='Image Y (um)', aspect='equal',
            title='Ch{} | mean {:.4g} W'.format(i, report['metrics']['mean_image_power_w'][i-1]))
    fig.colorbar(image, ax=list(axes), label='Mean image flux (W / um2)')
    fig.suptitle('Whole-beam first-order haze | independent objectives | FFT images are periodic')
    fig.savefig(output/'whole_haze_images.png', dpi=140)
    plt.close(fig)
    plot_statistics(output/'whole_haze_statistics.png', report['statistics'])
