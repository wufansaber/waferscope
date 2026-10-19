"""Direction annotations shared by upper-hemisphere far-field plots."""
import numpy as np


def mark_illumination_source(ax, incidence_deg):
    """Mark the sample-to-source direction for illumination arriving from -KY."""
    azimuth = -np.pi / 2
    radius = np.sin(np.deg2rad(incidence_deg))
    ax.scatter([azimuth], [radius], s=135, marker='*', facecolor='white',
               edgecolor='black', linewidth=1.2, zorder=8)
    ax.annotate('Laser from -KY\n{:g} deg incidence'.format(incidence_deg),
                xy=(azimuth, radius), xytext=(np.deg2rad(235), .55),
                ha='center', va='center', fontsize=8, color='black',
                bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='#333333', alpha=.9),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.2), zorder=9)


def mark_illumination_source_xy(ax, incidence_deg):
    """Mark source at global -KY in the SCATMECH Cartesian projection."""
    y = -np.sin(np.deg2rad(incidence_deg))
    ax.scatter([0], [y], s=105, marker='*', facecolor='white',
               edgecolor='black', linewidth=1.1, zorder=8)
    ax.annotate('Laser from -KY\n{:g} deg incidence'.format(incidence_deg),
                xy=(0, y), xytext=(-.7, -.72), ha='center', va='center',
                fontsize=8, color='black',
                bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='#333333', alpha=.9),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.2), zorder=9)
