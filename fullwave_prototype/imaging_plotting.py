"""Near field, angular far field and three separate image planes."""
import os
from pathlib import Path

import numpy as np

from .angular_plotting import angular_disk_edges, format_angular_disk
from .plot_directions import mark_illumination_source_xy


def _pyplot():
    os.environ.setdefault('MPLCONFIGDIR', str(Path('build/fem_validation/.matplotlib').resolve()))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _near(ax, data, report, fig):
    from .near_planes import saved_planes
    from .near_plotting import draw_plane, shared_norm
    planes = saved_planes(data)
    plane = next(p for p in planes if p.name == 'yz')
    im = draw_plane(ax, data, report, plane, shared_norm(data, planes))
    ax.set_title('1  NEAR FIELD | Total E, YZ (X = 0 nm)')
    fig.colorbar(im, ax=ax, label='|E / E_inc|^2 (sqrt colors)', shrink=.9)


def _far(ax, data, settings, fig, incidence_deg):
    from matplotlib.colors import PowerNorm
    theta, phi = data['theta'], data['phi']
    x_edges, y_edges = angular_disk_edges(theta, phi)
    im = ax.pcolormesh(x_edges, y_edges, data['dsc_total'],
                       norm=PowerNorm(.5, vmin=0), cmap='viridis', shading='flat')
    for index, (na, tilt) in enumerate(zip(settings['numerical_apertures'], settings['tilt_degrees']), 1):
        if na == 0:
            continue
        angle = np.linspace(0, 2*np.pi, 361)
        t = np.deg2rad(tilt)
        axis = np.array((np.sin(t), 0, np.cos(t)))
        u, v = np.array((np.cos(t), 0, -np.sin(t))), np.array((0, 1, 0))
        ring = np.sqrt(1-na*na)*axis+na*(np.cos(angle)[:, None]*u+np.sin(angle)[:, None]*v)
        ring[ring[:, 2] < 0, :2] = np.nan
        ax.plot(ring[:, 0], ring[:, 1], color=('#202020', '#ff7b54', '#efcb49')[index-1],
                linewidth=1.6, label='Ch{}: {:+g} deg, NA {:g}'.format(index, tilt, na))
    format_angular_disk(ax)
    ax.set_title('2  FAR FIELD | Upper hemisphere DSC')
    mark_illumination_source_xy(ax, incidence_deg)
    ax.legend(loc='upper left', bbox_to_anchor=(1.08, .3), fontsize=8)
    fig.colorbar(im, ax=ax, label='DSC (um^2 / sr)', pad=.12, shrink=.9)


def _image(ax, arrays, index, key, limit, norm, cmap):
    axis = arrays['ch{}_image_axis_um'.format(index)]
    value = arrays['ch{}_image_{}_w_um2'.format(index, key)]
    selected = abs(axis) <= limit+2*abs(axis[1]-axis[0])
    im = ax.pcolormesh(axis[selected], axis[selected], value[np.ix_(selected, selected)],
                       norm=norm, cmap=cmap, shading='auto')
    ax.set(xlim=(-limit, limit), ylim=(-limit, limit), xlabel="Image x' (um)", ylabel="Image y' (um)", aspect='equal')
    return im


def plot_pipeline(output, data, arrays, scattering, report):
    from matplotlib.colors import PowerNorm, TwoSlopeNorm
    plt = _pyplot()
    settings = report['settings']
    limit = 15*settings['magnification']/20
    vmax = max(arrays['ch{}_image_with_w_um2'.format(i)].max() for i in range(1, 4))
    norm = PowerNorm(.5, vmin=0, vmax=max(vmax, 1e-30))
    fig = plt.figure(figsize=(16, 10), layout='constrained')
    gs = fig.add_gridspec(2, 6)
    _near(fig.add_subplot(gs[0, :3]), data, scattering, fig)
    _far(fig.add_subplot(gs[0, 3:]), data, settings, fig,
         scattering['parameters']['incidence_deg'])
    image_axes = []
    for index, info in enumerate(report['channels'], 1):
        ax = fig.add_subplot(gs[1, (index-1)*2:index*2])
        image_axes.append(ax)
        im = _image(ax, arrays, index, 'with', limit, norm, 'magma')
        ax.set_title('3  IMAGE Ch{} | {:+g} deg, NA {:g}\nFull-plane P = {:.4g} pW'.format(
            index, info['tilt_deg'], info['na'], info['with']['image_power_w']*1e12), fontsize=11)
    fig.colorbar(im, ax=image_axes, label='I (W / um^2), shared absolute scale', shrink=.9)
    p = scattering['parameters']
    fig.suptitle('FULL-WAVE DARK-FIELD PIPELINE | {} nm, particle {} nm | {}x, f_obj={} mm, f_tube={} mm\n'
                 'Experimental FEM + configured optics; independent channels; far field / images exclude specular beam'.format(
                     p['wavelength_nm'], p['diameter_nm'], settings['magnification'],
                     settings['focal_length_mm'], report['tube_focal_length_mm']), fontsize=12)
    fig.savefig(Path(output)/'near_far_image.png', dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(3, 4, figsize=(17, 12), layout='constrained')
    for column, (key, label) in enumerate((('with', 'With particle'), ('without', 'Roughness background'),
                                         ('change', 'Changed-field flux'), ('interference', 'Signed interference'))):
        maximum = max(abs(arrays['ch{}_image_{}_w_um2'.format(i, key)]).max() for i in range(1, 4))
        maximum = max(maximum, 1e-30)
        cnorm = TwoSlopeNorm(vmin=-maximum, vcenter=0, vmax=maximum) if key == 'interference' else PowerNorm(.5, vmin=0, vmax=maximum)
        for row in range(3):
            im = _image(axes[row, column], arrays, row+1, key, limit, cnorm,
                        'RdBu_r' if key == 'interference' else 'magma')
            axes[row, column].set_title('Ch{} | {}'.format(row+1, label), fontsize=10)
        fig.colorbar(im, ax=axes[:, column].tolist(), label='W / um^2', shrink=.7, location='bottom')
    fig.suptitle('Coherent image decomposition | Total = background + changed-field flux + signed interference\n'
                 'Each column shares an absolute scale across channels; columns have separate labeled scales', fontsize=12)
    fig.savefig(Path(output)/'image_decomposition.png', dpi=150)
    plt.close(fig)
    plot_stages(output, data, arrays, scattering, report)


def plot_stages(output, data, arrays, scattering, report):
    from matplotlib.colors import PowerNorm
    plt = _pyplot()
    settings = report['settings']
    from .near_plotting import plot_near_planes
    plot_near_planes(output, data, scattering)
    fig = plt.figure(figsize=(10, 7), layout='constrained')
    _far(fig.add_subplot(), data, settings, fig,
         scattering['parameters']['incidence_deg'])
    fig.savefig(Path(output)/'far_field.png', dpi=150)
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), layout='constrained')
    maximum = max(arrays['ch{}_image_with_w_um2'.format(i)].max() for i in range(1, 4))
    norm = PowerNorm(.5, vmin=0, vmax=max(maximum, 1e-30))
    for i, (ax, channel) in enumerate(zip(axes, report['channels']), 1):
        im = _image(ax, arrays, i, 'with', 15*settings['magnification']/20, norm, 'magma')
        ax.set_title('Ch{} | {:+g} deg | NA {:g}\n{:.4g} pW'.format(
            i, channel['tilt_deg'], channel['na'], channel['with']['image_power_w']*1e12))
    fig.colorbar(im, ax=axes.tolist(), label='W / um^2; shared absolute scale', shrink=.8)
    fig.savefig(Path(output)/'image_planes.png', dpi=150)
    plt.close(fig)
    plot_image_channels(output, arrays, report, norm)


def plot_image_channels(output, arrays, report, norm):
    """Export independent channels without renormalizing their relative powers."""
    plt = _pyplot()
    settings = report['settings']
    for index, channel in enumerate(report['channels'], 1):
        fig, ax = plt.subplots(figsize=(8, 7), layout='constrained')
        im = _image(ax, arrays, index, 'with', 15*settings['magnification']/20, norm, 'magma')
        ax.set_title('Ch{} | {:+g} deg | NA {:g} | {}x\nFull-plane power: {:.4g} pW'.format(
            index, channel['tilt_deg'], channel['na'], settings['magnification'],
            channel['with']['image_power_w']*1e12))
        fig.colorbar(im, ax=ax, label='W / um^2; shared scale (sqrt colors)', shrink=.85)
        fig.suptitle('IMAGE PLANE | Total scattered light with particle')
        fig.savefig(Path(output)/('image_channel_{}.png'.format(index)), dpi=150)
        plt.close(fig)
