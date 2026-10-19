"""Unnormalized rough-patch fields and coherent scattering decomposition."""
import os
from pathlib import Path

import numpy as np


def plot_pair(output, data, report):
    os.environ.setdefault('MPLCONFIGDIR', str(Path('build/fem_validation/.matplotlib').resolve()))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import PowerNorm
    from matplotlib.tri import Triangulation
    from .defects import draw_section, description
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), layout='constrained')
    triangles = data['surface_triangles_um']*1000
    triangulation = Triangulation(triangles[..., 0].ravel(), triangles[..., 1].ravel(),
                                  np.arange(triangles.shape[0]*3).reshape(-1, 3))
    limit = max(abs(triangles[..., 2]).max(), 1e-10)
    im = axes[0, 0].tripcolor(triangulation, triangles[..., 2].ravel(), shading='gouraud',
                             cmap='RdBu_r', vmin=-limit, vmax=limit)
    fig.colorbar(im, ax=axes[0, 0], label='Height (nm)')
    axes[0, 0].set(title='Finite rough patch', xlabel='X (nm)', ylabel='Y (nm)', aspect='equal')
    p = report['parameters']
    for ax, key, title in ((axes[0, 1], 'near_e_with', 'Total field with particle'),
                           (axes[0, 2], 'near_e_change', 'Field change caused by particle')):
        value = np.sum(abs(data[key])**2, axis=-1)
        im = ax.pcolormesh(data['near_y_um']*1000, data['near_z_um']*1000, value,
                           norm=PowerNorm(.5, vmin=0), cmap='magma', shading='auto')
        fig.colorbar(im, ax=ax, label='|E / E_inc|^2 (sqrt color scale)')
        surface = data['surface_vertices_um']
        # Default odd grid contains the x=0 row; interpolate for even grids.
        y = surface[0, :, 1]*1000
        z = [np.interp(0, surface[:, 0, 0], surface[:, j, 2])*1000 for j in range(len(y))]
        ax.plot(y, z, color='cyan', linewidth=.8)
        draw_section(ax, p, 'yz', linewidth=.8)
        ax.set(title=title, xlabel='Y (nm)', ylabel='Z (nm)', aspect='equal')
    colors = dict(total='#174c78', background='#777777', change_field='#c45714', interference='#16816e')
    labels = dict(total='With particle', background='Without particle', change_field='Field-change power',
                  interference='Interference (signed)')
    for ax, angle in ((axes[1, 0], 90), (axes[1, 1], 270)):
        index = np.argmin(abs(np.rad2deg(data['phi'])-angle))
        for key in colors:
            ax.plot(np.rad2deg(data['theta']), data['dsc_'+key][:, index], color=colors[key],
                    label=labels[key], linewidth=1.5)
        ax.axhline(0, color='black', linewidth=.5)
        ax.set(title='Far field: phi = {} deg'.format(angle), xlabel='Theta from +Z (deg)',
               ylabel='DSC (um^2 / sr)', xlim=(0, 90))
        ax.grid(alpha=.2)
    axes[1, 0].legend(fontsize=8)
    ax = axes[1, 2]
    integrals = report['metrics']['hemisphere_integrals_um2']
    ax.bar(range(4), [integrals[k] for k in colors], color=list(colors.values()))
    ax.set_xticks(range(4), ['With\nparticle', 'Without\nparticle', 'Field\nchange', 'Interference'])
    ax.axhline(0, color='black', linewidth=.5)
    ax.set(title='Upper-hemisphere integral', ylabel='Integrated DSC (um^2)')
    sample_rms = report['surface_statistics']['triangle_rms_nm']
    fig.suptitle('EXPERIMENTAL | {} nm, {}, rough sample RMS {:.3g} nm, defect bottom z={} nm\n'
                 'Same geometry and phase reference; total = background + change-field power + interference'.format(
                     p['wavelength_nm'], description(p), sample_rms, p['gap_nm']), fontsize=12)
    fig.savefig(Path(output)/'rough_comparison.png', dpi=150)
    plt.close(fig)
