"""Comparable orthogonal near-field plots from phase-resolved FEM samples."""
from pathlib import Path

import numpy as np

from .near_planes import PLANE_AXES, electric_field, saved_planes


def surface_segments(triangles, fixed_axis, offset):
    """Intersect actual mesh triangles, including diagonals on even rough grids."""
    segments = []
    for triangle in triangles:
        distances = triangle[:, fixed_axis]-offset
        hits = []
        for i, j in ((0, 1), (1, 2), (2, 0)):
            if abs(distances[i]) < 1e-12:
                hits.append(triangle[i])
            if distances[i]*distances[j] < 0:
                fraction = distances[i]/(distances[i]-distances[j])
                hits.append(triangle[i]+fraction*(triangle[j]-triangle[i]))
        if hits:
            hits = np.unique(np.round(hits, 12), axis=0)
            if len(hits) == 2:
                segments.append(hits)
    return segments


def draw_plane(ax, data, report, plane, norm):
    from matplotlib.collections import LineCollection
    from matplotlib.patches import Circle
    horizontal, vertical, fixed = PLANE_AXES[plane.name]
    value = np.sum(abs(electric_field(data, plane))**2, axis=-1)
    im = ax.pcolormesh(plane.u*1000, plane.v*1000, value,
                       norm=norm, cmap='magma', shading='auto')
    p = report['parameters']
    center = np.array((0., 0., p['diameter_nm']/2+p['gap_nm']))
    radius_squared = (p['diameter_nm']/2)**2-(plane.offset*1000-center[fixed])**2
    if radius_squared > 0:
        ax.add_patch(Circle(center[[horizontal, vertical]], np.sqrt(radius_squared),
                            fill=False, edgecolor='cyan', linewidth=1))
    if fixed != 2:
        ax.axhline(0, color='cyan', linewidth=.6, linestyle=':')
    segments = surface_segments(data['surface_triangles_um'], fixed, plane.offset)
    if segments:
        projected = [segment[:, [horizontal, vertical]]*1000 for segment in segments]
        ax.add_collection(LineCollection(projected, colors='cyan', linewidths=1))
    ax.set(title='{} | {} = {:g} nm'.format(plane.name.upper(), 'XYZ'[fixed], plane.offset*1000),
           xlabel='{} (nm)'.format('XYZ'[horizontal]), ylabel='{} (nm)'.format('XYZ'[vertical]),
           aspect='equal', xlim=(plane.u[0]*1000, plane.u[-1]*1000),
           ylim=(plane.v[0]*1000, plane.v[-1]*1000))
    return im


def shared_norm(data, planes):
    from matplotlib.colors import PowerNorm
    maximum = max(float(np.max(np.sum(abs(electric_field(data, p))**2, axis=-1))) for p in planes)
    return PowerNorm(.5, vmin=0, vmax=max(maximum, 1e-30))


def plot_near_planes(output, data, report):
    from .imaging_plotting import _pyplot
    plt = _pyplot()
    planes = saved_planes(data)
    norm = shared_norm(data, planes)
    label = '|E / E_inc|^2; shared scale (sqrt colors)'
    fig, axes = plt.subplots(1, len(planes), figsize=(5.5*len(planes), 5.5),
                             squeeze=False, layout='constrained')
    for ax, plane in zip(axes[0], planes):
        im = draw_plane(ax, data, report, plane, norm)
    fig.colorbar(im, ax=axes.ravel().tolist(), label=label, shrink=.8)
    fig.suptitle('NEAR FIELD | Total electric field with particle')
    fig.savefig(Path(output)/'near_field.png', dpi=150)
    plt.close(fig)
    for plane in planes:
        fig, ax = plt.subplots(figsize=(8, 7), layout='constrained')
        im = draw_plane(ax, data, report, plane, norm)
        fig.colorbar(im, ax=ax, label=label, shrink=.9)
        fig.suptitle('NEAR FIELD | Total electric field with particle')
        fig.savefig(Path(output)/('near_field_'+plane.name+'.png'), dpi=150)
        plt.close(fig)
