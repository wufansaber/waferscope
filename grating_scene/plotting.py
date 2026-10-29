"""Unnormalized periodic-scene figures with explicit physical coordinates."""
import numpy as np
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm


def positive_norm(arrays, logarithmic=False):
    values = np.concatenate([np.asarray(a).ravel() for a in arrays])
    upper = float(np.max(values)) if values.size else 0.
    positive = values[values > 0]
    if logarithmic and len(positive) and upper > 0:
        lower = max(float(positive.min()), upper*1e-8)
        return LogNorm(lower, max(upper, lower*1.01))
    return Normalize(0., upper if upper > 0 else 1.)


def geometry_figure(fig, config):
    fig.clear()
    ax = fig.add_subplot(111, projection='3d')
    s = config.scene
    x = np.linspace(-s.lx_um/2, s.lx_um/2, 61)
    y = np.linspace(-s.ly_um/2, s.ly_um/2, 121)
    xx, yy = np.meshgrid(x, y)
    offset = (yy+s.pitch_um/2) % s.pitch_um-s.pitch_um/2
    zz = np.where(abs(offset) <= s.ridge_width_um/2, 0., -s.depth_um)
    if s.kind == 'flat':
        zz[:] = 0.
    ax.plot_surface(xx, yy, zz, color='#5D829A', linewidth=0, alpha=.8)
    if config.particle:
        phi, theta = np.meshgrid(np.linspace(0, 2*np.pi, 33), np.linspace(0, np.pi, 25))
        r = s.diameter_um/2
        ax.plot_surface(r*np.sin(theta)*np.cos(phi), r*np.sin(theta)*np.sin(phi),
                        s.gap_um+r+r*np.cos(theta), color='#C96845', linewidth=0)
    angle = np.deg2rad(s.incidence_deg)
    incident = np.array([0., np.sin(angle), -np.cos(angle)])
    end = np.array([0., 0., s.gap_um+s.diameter_um/2])
    length = .45*max(s.lx_um, s.ly_um)
    start = end-length*incident
    ax.quiver(*start, *(length*incident), color='#274E8B', linewidth=1.5, arrow_length_ratio=.15)
    ax.set(xlabel='x (um)', ylabel='y (um)', zlabel='z (um)',
           title='Periodic cell: source from -KY, propagation toward +KY')
    ax.set_box_aspect((s.lx_um, s.ly_um, max(s.depth_um+s.cap_top_um, .2)))
    fig.tight_layout()


def near_figure(fig, near, plane, logarithmic=False):
    fig.clear()
    points = near[plane+'_points_um']
    axes_indices = {'xy': (0, 1), 'xz': (0, 2), 'yz': (1, 2)}[plane]
    a, b = axes_indices
    extent = [points[0, 0, a], points[0, -1, a], points[0, 0, b], points[-1, 0, b]]
    values = [np.sum(abs(near[plane+'_'+label+'_electric'])**2, axis=-1)
              for label in ('reference', 'sample', 'delta')]
    shared = positive_norm(values[:2], logarithmic)
    for i, (label, data) in enumerate(zip(('Reference: total E', 'Sample: total E', 'Difference field'), values), 1):
        ax = fig.add_subplot(1, 3, i)
        im = ax.imshow(data, origin='lower', extent=extent, aspect='auto', cmap='viridis',
                       norm=shared if i < 3 else positive_norm([data], logarithmic), interpolation='nearest')
        ax.set(title=label, xlabel='xyz'[a]+' (um)', ylabel='xyz'[b]+' (um)')
        fig.colorbar(im, ax=ax, shrink=.75, label='|E/E0|^2' if i < 3 else '|(Es-Er)/E0|^2')
    fig.suptitle(plane.upper()+' | native total fields; separate difference-field scale')
    fig.tight_layout()


def far_figure(fig, far, config, logarithmic=False, three_dimensional=False):
    fig.clear()
    valid = far['propagating']
    rays = far['directions'][valid]
    values = [far['efficiency_'+label][valid] for label in ('reference', 'sample', 'delta')]
    norm = positive_norm(values, logarithmic)
    for i, (label, power) in enumerate(zip(('Reference', 'Sample', 'Difference field'), values), 1):
        ax = fig.add_subplot(1, 3, i, projection='3d' if three_dimensional else None)
        if three_dimensional:
            im = ax.scatter(*rays.T, c=power, s=30, cmap='viridis', norm=norm)
            for na, tilt in zip(config.numerical_apertures, config.tilts_deg):
                if na == 0:
                    continue
                t, p = np.deg2rad(tilt), np.linspace(0, 2*np.pi, 120)
                axis, u, v = np.array([np.sin(t), 0, np.cos(t)]), np.array([np.cos(t), 0, -np.sin(t)]), np.array([0, 1, 0])
                ring = np.sqrt(1-na*na)*axis+na*(np.cos(p)[:, None]*u+np.sin(p)[:, None]*v)
                ring = ring[ring[:, 2] > 0]
                ax.plot(*ring.T, lw=1)
                for ray in ring[::max(1, len(ring)//8)]:
                    ax.plot([0, ray[0]], [0, ray[1]], [0, ray[2]], color='.6', lw=.4)
            ax.set(xlim=(-1, 1), ylim=(-1, 1), zlim=(0, 1), zlabel='kz/k')
            ax.set_box_aspect((2, 2, 1))
        else:
            im = ax.scatter(rays[:, 0], rays[:, 1], c=power, s=45, cmap='viridis', norm=norm)
            p = np.linspace(0, 2*np.pi, 181)
            ax.plot(np.cos(p), np.sin(p), color='.5', lw=.7)
            for na, tilt in zip(config.numerical_apertures, config.tilts_deg):
                if na:
                    t = np.deg2rad(tilt)
                    ring = np.column_stack((np.sqrt(1-na*na)*np.sin(t)+na*np.cos(p)*np.cos(t),
                                            na*np.sin(p), np.sqrt(1-na*na)*np.cos(t)-na*np.cos(p)*np.sin(t)))
                    ring[ring[:, 2] <= 0] = np.nan
                    ax.plot(ring[:, 0], ring[:, 1], lw=1)
            ax.set(xlim=(-1.05, 1.05), ylim=(-1.05, 1.05), aspect='equal')
        ax.set(title=label, xlabel='kx/k', ylabel='ky/k')
        fig.colorbar(im, ax=ax, shrink=.65, label='Order efficiency')
    fig.suptitle('Discrete periodic orders, not continuous DSC | source: -KY')
    fig.tight_layout()


def image_figure(fig, images, collection, logarithmic=False):
    fig.clear()
    x, y = images['image_x_um'], images['image_y_um']
    dx, dy = x[1]-x[0], y[1]-y[0]
    extent = [x[0]-dx/2, x[-1]+dx/2, y[0]-dy/2, y[-1]+dy/2]
    for i in range(1, 4):
        ref, sample, delta = [images['ch{}_{}_w_um2'.format(i, key)]
                              for key in ('reference', 'sample', 'difference')]
        shared = positive_norm([ref, sample], logarithmic)
        limit = max(float(abs(delta).max()), 1e-300)
        for j, (label, data) in enumerate(zip(('Reference', 'Sample', 'Sample - reference'), (ref, sample, delta)), 1):
            ax = fig.add_subplot(3, 3, (i-1)*3+j)
            norm = shared if j < 3 else TwoSlopeNorm(vmin=-limit, vcenter=0., vmax=limit)
            im = ax.imshow(data, origin='lower', extent=extent, aspect='auto',
                           norm=norm, cmap='viridis' if j < 3 else 'RdBu_r', interpolation='nearest')
            state = '' if collection['channels'][i-1]['enabled'] else ' OFF'
            ax.set(title='CH{}{} {}'.format(i, state, label), xlabel='Image x (um)', ylabel='Image y (um)')
            fig.colorbar(im, ax=ax, shrink=.7, label='W/um^2')
    fig.suptitle('Ideal wafer-conjugate periodic Jones relay; not camera data')
    fig.tight_layout()


def save_figures(config, near, far, images, collection, output):
    from matplotlib.figure import Figure
    renderers = [('geometry.png', (8, 5), lambda f: geometry_figure(f, config)),
        ('far_field_2d.png', (12, 4.2), lambda f: far_figure(f, far, config, True)),
        ('far_field_3d.png', (12, 4.5), lambda f: far_figure(f, far, config, True, True)),
        ('images.png', (13, 11), lambda f: image_figure(f, images, collection, True))]
    for plane in ('xy', 'xz', 'yz'):
        if plane+'_points_um' in near:
            renderers.append(('near_'+plane+'.png', (12, 4), lambda f, p=plane: near_figure(f, near, p, True)))
    for name, size, render in renderers:
        fig = Figure(figsize=size)
        render(fig)
        fig.savefig(output/name, dpi=140)
