"""Static benchmark figures with absolute DSC and explicit physical assumptions."""
from pathlib import Path
import os
import numpy as np


def plot_result(output, p, metrics, reference=None):
    cache = Path(__file__).resolve().parents[1]/'build'/'fem_validation'/'.matplotlib'
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR', str(cache))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import PowerNorm
    from matplotlib.patches import Circle
    output = Path(output)
    with np.load(output/'fields.npz') as data:
        e = data['near_e']
        esc = data['near_scattered_e']
        y, z = data['near_y_um'], data['near_z_um']
        theta, phi = data['theta'], data['phi']
        dsc = data['dsc_um2_sr']
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout='constrained')
    for ax, values, title in zip(axes[0], (e, esc), ('Total field |E/Einc|^2', 'Scattered field |Esc/Einc|^2')):
        intensity = np.sum(abs(values)**2, axis=-1)
        im = ax.imshow(intensity, origin='lower', extent=(y[0]*1000, y[-1]*1000, z[0]*1000, z[-1]*1000),
                       cmap='inferno', norm=PowerNorm(.5), aspect='equal')
        ax.axhline(0, color='cyan', lw=.7)
        ax.add_patch(Circle((0, p.diameter_nm/2+p.gap_nm), p.diameter_nm/2,
                            fill=False, edgecolor='cyan', linewidth=.8))
        ax.set(title=title, xlabel='Y (nm)', ylabel='Z (nm)')
        fig.colorbar(im, ax=ax)
    for ax, azimuth in zip(axes[1], (0, 90)):
        index = int(np.argmin(abs(np.rad2deg(phi)-azimuth)))
        ax.plot(np.rad2deg(theta), dsc[:, index], label='FEM')
        if reference is not None:
            ax.plot(np.rad2deg(theta), reference[:, index], '--', label='SCATMECH lmax=14')
        ax.set(title='Far field, azimuth {} deg'.format(azimuth), xlabel='Polar angle from +Z (deg)',
               ylabel='DSC (um^2/sr)')
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        ax.grid(alpha=.25)
        ax.legend()
    title = 'Experimental FEM | {:.0f} nm laser | {:.0f} nm sphere | gap {:.3g} nm | {} polarization'.format(
        p.wavelength_nm, p.diameter_nm, p.gap_nm, p.polarization.upper())
    if reference is not None:
        title += '\nHemisphere power discrepancy {:.2%} (not a full error estimate)'.format(metrics['hemisphere_relative_error'])
    fig.suptitle(title, fontsize=12)
    fig.savefig(output/'comparison.png', dpi=160)
    plt.close(fig)
