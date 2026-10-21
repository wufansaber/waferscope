"""Geometry and absolute field plots for the experimental finite 3D grating."""
import numpy as np


def geometry_faces(g):
    from .finite_grating import trench_vertices
    faces, kind = [], []
    def append(points, category):
        faces.append(np.asarray(points))
        kind.append(category)
    h, length = g.half_width_um, g.length_um/2
    edges = [-h]+[v for c in g.centers for v in (c-g.top_width_um/2, c+g.top_width_um/2)]+[h]
    for a, b in zip(edges[:-1], edges[1:]):
        inside = any(a < c < b for c in g.centers)
        if not inside:
            append([[-length, a, 0], [length, a, 0], [length, b, 0], [-length, b, 0]], 'top')
        for lo, hi in ((-h, -length), (length, h)):
            append([[lo, a, 0], [hi, a, 0], [hi, b, 0], [lo, b, 0]], 'top')
    for center in g.centers:
        vertices = trench_vertices(g, center)
        for ids in ([4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]):
            append(vertices[ids], 'groove')
    return faces, kind


def plot_geometry(output, g):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    faces, kinds = geometry_faces(g)
    fig = plt.figure(figsize=(9, 6), constrained_layout=True)
    ax = fig.add_subplot(111, projection='3d')
    colors = ['#b4c5ce' if k == 'top' else '#677e8c' for k in kinds]
    ax.add_collection3d(Poly3DCollection(faces, facecolors=colors, edgecolors='#364952', linewidths=.5))
    h = g.half_width_um
    ax.set(xlim=(-h, h), ylim=(-h, h), zlim=(-g.depth_um*1.15, g.depth_um*.3),
           xlabel='Line x (um)', ylabel='Transverse y (um)', zlabel='Height z (um)')
    ax.set_box_aspect((2*h, 2*h, g.depth_um*2.5))
    ax.view_init(elev=33, azim=-58)
    ax.set_title('{} finite trenches | pitch {:g} nm | length {:g} um\n'
        'Top / bottom width {:g} / {:g} nm | depth {:g} nm'.format(
        g.trench_count, g.pitch_um*1000, g.length_um, g.top_width_um*1000,
        g.bottom_width_um*1000, g.depth_um*1000))
    fig.savefig(output/'grating_geometry_3d.png', dpi=180)
    plt.close(fig)


def plot_fields(output, data, report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    values = [data['near_'+name+'_intensity'] for name in ('xy', 'xz', 'yz')]
    vmax = max(float(v.max()) for v in values)
    norm = LogNorm(max(vmax*1e-5, 1e-10), vmax)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.6), constrained_layout=True)
    for ax, name, value in zip(axes, ('xy', 'xz', 'yz'), values):
        artist = ax.pcolormesh(data['near_'+name+'_u'], data['near_'+name+'_v'], value,
                              norm=norm, cmap='inferno', shading='auto')
        fixed = next(axis for axis in 'xyz' if axis not in name)
        ax.set(title='{} | {} = {:g} um'.format(name.upper(), fixed, float(data['near_'+name+'_offset_um'])),
               xlabel=name[0]+' (um)', ylabel=name[1]+' (um)', aspect='equal')
    fig.colorbar(artist, ax=list(axes), shrink=.7, label='Total |E / E_inc|^2 (shared log scale)')
    fig.suptitle('Finite 3D grating | total near field | numerical verification pending', fontsize=12)
    fig.savefig(output/'near_field_xy_xz_yz.png', dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
    artist = ax.pcolormesh(data['disk_axis'], data['disk_axis'], data['disk_dsc'],
                          cmap='viridis', shading='auto', vmin=0)
    angle = np.linspace(0, 2*np.pi, 361)
    ax.plot(np.cos(angle), np.sin(angle), color='#777777', lw=.8)
    incident_y = -np.sin(np.deg2rad(report['parameters']['incidence_deg']))
    ax.plot(0, incident_y, marker='*', color='white', markeredgecolor='black', ms=10)
    ax.annotate('Source on -KY side', xy=(0, incident_y), xytext=(-.85, -.72), fontsize=8,
                bbox=dict(facecolor='white', edgecolor='#aaaaaa', alpha=.9),
                arrowprops=dict(arrowstyle='->', color='black'))
    ax.set(xlabel='Global kx/k', ylabel='Global ky/k', aspect='equal',
           title='Finite grating scattered DSC\nUnvalidated FEM result; no specular background')
    fig.colorbar(artist, ax=ax, label='DSC (um^2 / sr)')
    fig.savefig(output/'far_field_grating.png', dpi=180)
    plt.close(fig)
