"""Near field, angular far field and three separate image planes."""
import os
from pathlib import Path

import numpy as np


def _pyplot():
    os.environ.setdefault('MPLCONFIGDIR', str(Path('build/fem_validation/.matplotlib').resolve()))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _near(ax, data, report, fig):
    from matplotlib.colors import PowerNorm
    from matplotlib.patches import Circle
    value = np.sum(abs(data['near_e_with'])**2, axis=-1)
    im = ax.pcolormesh(data['near_y_um']*1000, data['near_z_um']*1000, value,
                       norm=PowerNorm(.5, vmin=0), cmap='magma', shading='auto')
    p = report['parameters']
    ax.add_patch(Circle((0, p['diameter_nm']/2+p['gap_nm']), p['diameter_nm']/2,
                        fill=False, edgecolor='cyan', linewidth=1))
    vertices = data['surface_vertices_um']
    j = np.argmin(abs(vertices[:, 0, 0]))
    ax.plot(vertices[j, :, 1]*1000, vertices[j, :, 2]*1000, color='cyan', linewidth=1)
    ax.set(title='1  NEAR FIELD | Total E, YZ section', xlabel='Y (nm)', ylabel='Z (nm)', aspect='equal')
    fig.colorbar(im, ax=ax, label='|E / E_inc|^2', shrink=.9)


def _far(ax, data, settings, fig):
    from matplotlib.colors import PowerNorm
    theta, phi = data['theta'], data['phi']
    theta_edges = np.r_[0, (theta[:-1]+theta[1:])/2, np.pi/2]
    dphi = phi[1]-phi[0]
    phi_edges = np.r_[phi-dphi/2, phi[-1]+dphi/2]
    im = ax.pcolormesh(phi_edges, np.sin(theta_edges), data['dsc_total'],
                       norm=PowerNorm(.5, vmin=0), cmap='viridis', shading='flat')
    for index, (na, tilt) in enumerate(zip(settings['numerical_apertures'], settings['tilt_degrees']), 1):
        if na == 0:
            continue
        angle = np.linspace(0, 2*np.pi, 361)
        t = np.deg2rad(tilt)
        axis = np.array((np.sin(t), 0, np.cos(t)))
        u, v = np.array((np.cos(t), 0, -np.sin(t))), np.array((0, 1, 0))
        ring = np.sqrt(1-na*na)*axis+na*(np.cos(angle)[:, None]*u+np.sin(angle)[:, None]*v)
        azimuth, radius = np.arctan2(ring[:, 1], ring[:, 0]), np.hypot(ring[:, 0], ring[:, 1])
        radius[ring[:, 2] < 0] = np.nan
        ax.plot(azimuth, radius, color=('#202020', '#ff7b54', '#efcb49')[index-1],
                linewidth=1.6, label='Ch{}: {:+g} deg, NA {:g}'.format(index, tilt, na))
    ax.set(ylim=(0, 1), title='2  FAR FIELD | Upper hemisphere DSC')
    ax.set_yticks([.5, 1], labels=['0.5', '1'])
    ax.set_rlabel_position(70)
    ax.grid(alpha=.3)
    ax.legend(loc='upper left', bbox_to_anchor=(1.13, .3), fontsize=8)
    fig.colorbar(im, ax=ax, label='um^2 / sr; radial coordinate = sin(theta)', pad=.16, shrink=.9)


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
    _far(fig.add_subplot(gs[0, 3:], projection='polar'), data, settings, fig)
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
                 'Experimental FEM + ideal optics; independent channels; far field / images exclude specular beam'.format(
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
    fig, ax = plt.subplots(figsize=(9, 7), layout='constrained')
    _near(ax, data, scattering, fig)
    fig.savefig(Path(output)/'near_field.png', dpi=150)
    plt.close(fig)
    fig = plt.figure(figsize=(10, 7), layout='constrained')
    _far(fig.add_subplot(projection='polar'), data, settings, fig)
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
