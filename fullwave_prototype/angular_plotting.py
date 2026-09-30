"""Global kx/k, ky/k projection for saved upper-hemisphere angular grids."""

import numpy as np
from matplotlib.patches import Circle


def angular_disk_edges(theta, phi):
    """Map angular cell edges to the same global XY disk as SCATMECH."""
    theta = np.asarray(theta)
    phi = np.asarray(phi)
    theta_edges = np.r_[0, (theta[:-1]+theta[1:])/2, np.pi/2]
    dphi = phi[1]-phi[0]
    phi_edges = np.r_[phi-dphi/2, phi[-1]+dphi/2]
    polar, azimuth = np.meshgrid(theta_edges, phi_edges, indexing='ij')
    return np.sin(polar)*np.cos(azimuth), np.sin(polar)*np.sin(azimuth)


def format_angular_disk(ax):
    """Use SCATMECH's normalized global transverse-wavevector axes."""
    ax.add_patch(Circle((0, 0), 1, fill=False, color='#202020', lw=1.2, zorder=7))
    ax.set(xlim=(-1.05, 1.05), ylim=(-1.05, 1.05), aspect='equal',
           xlabel=r'$k_x/k$', ylabel=r'$k_y/k$')
    ax.set_xticks((-1, -.5, 0, .5, 1))
    ax.set_yticks((-1, -.5, 0, .5, 1))
    ax.grid(alpha=.18)
