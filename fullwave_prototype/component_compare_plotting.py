"""Far-field and image figures for independently scaled FEM components."""
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


def _pupil_rings(ax, settings):
    colors = ('#f4d35e', '#ee964b', '#f95738')
    for index, (na, tilt) in enumerate(zip(settings['numerical_apertures'],
                                           settings['tilt_degrees'])):
        if na == 0:
            continue
        angle = np.linspace(0, 2*np.pi, 361)
        t = np.deg2rad(tilt)
        axis = np.array((np.sin(t), 0, np.cos(t)))
        u = np.array((np.cos(t), 0, -np.sin(t)))
        v = np.array((0, 1, 0))
        ring = np.sqrt(1-na*na)*axis+na*(np.cos(angle)[:, None]*u+np.sin(angle)[:, None]*v)
        ring[ring[:, 2] < 0, :2] = np.nan
        ax.plot(ring[:, 0], ring[:, 1],
                color=colors[index], linewidth=1.2, label=f'Ch{index+1}')


def _component_norm(values, logarithmic):
    from matplotlib.colors import LogNorm, PowerNorm
    maximum = max(float(np.max(values)), 1e-30)
    return (LogNorm(vmin=maximum*1e-8, vmax=maximum) if logarithmic else
            PowerNorm(.5, vmin=0, vmax=maximum))


def _far(output, far, report, logarithmic=False):
    plt = _pyplot()
    theta, phi = far['theta'], far['phi']
    x_edges, y_edges = angular_disk_edges(theta, phi)
    fig = plt.figure(figsize=(17, 5.8), layout='constrained')
    gs = fig.add_gridspec(1, 3, width_ratios=(1, 1, .8))
    axes = [fig.add_subplot(gs[0, index]) for index in (0, 1)]
    totals = report['far_hemisphere']
    for ax, key, label in zip(axes, ('particle', 'haze'), ('Particle only', 'Haze only')):
        im = ax.pcolormesh(x_edges, y_edges, far[key+'_dsc'],
                           norm=_component_norm(far[key+'_dsc'], logarithmic),
                           cmap='viridis', shading='flat', rasterized=True)
        _pupil_rings(ax, report['imaging_settings'])
        format_angular_disk(ax)
        mark_illumination_source_xy(ax, report['parameters']['incidence_deg'])
        ax.set_title(f'{label}\nHemisphere DSC = {totals[key+"_dsc_um2"]:.4g} um^2')
        fig.colorbar(im, ax=ax, label='DSC (um^2 / sr)', shrink=.75, pad=.08)
    ax = fig.add_subplot(gs[0, 2])
    powers = [totals[key+'_equivalent_power_w']*1e12 for key in ('particle', 'haze')]
    ax.bar(('Particle', 'Haze'), powers, color=('#187b87', '#d27831'))
    for index, power in enumerate(powers):
        ax.text(index, power, f'{power:.4g} pW', ha='center', va='bottom')
    ax.set(ylabel='Upper-hemisphere equivalent power (pW)', title='Hemisphere integral')
    ax.grid(axis='y', alpha=.2)
    fig.suptitle('Independent FEM experiments | global kx/k, ky/k | separate color scales'+
                 (' | logarithmic colors' if logarithmic else ''))
    fig.savefig(output/('far_components_log.png' if logarithmic else 'far_components.png'), dpi=160)
    plt.close(fig)


def _images(output, fields, report, logarithmic):
    plt = _pyplot()
    active = [channel for channel in report['channels'] if channel['na'] > 0]
    if not active:
        active = report['channels']
    norms = {kind: _component_norm([fields[f'ch{channel["channel"]}_{kind}_w_um2'].max()
                                   for channel in active], logarithmic)
             for kind in ('particle', 'haze')}
    fig, axes = plt.subplots(len(active), 2, figsize=(12, 4*len(active)),
                             sharex=True, sharey=True, squeeze=False, layout='constrained')
    limit = 15*report['imaging_settings']['magnification']/20
    images = {}
    for row, channel in enumerate(active):
        index = channel['channel']
        axis = fields[f'ch{index}_axis_um']
        selected = abs(axis) <= limit+2*abs(axis[1]-axis[0])
        for column, kind in enumerate(('particle', 'haze')):
            ax = axes[row, column]
            power = channel[kind+'_power_w']*1e12
            value = fields[f'ch{index}_{kind}_w_um2']
            images[kind] = ax.pcolormesh(axis[selected], axis[selected],
                                         value[np.ix_(selected, selected)], norm=norms[kind],
                                         cmap='magma', shading='auto', rasterized=True)
            ax.set(xlim=(-limit, limit), ylim=(-limit, limit), aspect='equal',
                   title=f'Ch{index} {kind} | {channel["tilt_deg"]:+g} deg, NA {channel["na"]:g}\n'
                         f'Full-image P = {power:.4g} pW')
        axes[row, 0].set_ylabel("Image y' (um)")
    for ax in axes[-1]:
        ax.set_xlabel("Image x' (um)")
    for column, kind in enumerate(('particle', 'haze')):
        fig.colorbar(images[kind], ax=axes[:, column].tolist(),
                     label=f'{kind.capitalize()} I (W / um^2)', shrink=.7)
    fig.suptitle('Particle on smooth wafer vs rough wafer without particle'+
                 (' | logarithmic colors' if logarithmic else ' | square-root colors')+
                 ' | separate component scales')
    fig.savefig(output/('image_components_log.png' if logarithmic else 'image_components.png'), dpi=160)
    plt.close(fig)


def plot_comparison(output, far, fields, report):
    output = Path(output)
    _far(output, far, report)
    _far(output, far, report, logarithmic=True)
    _images(output, fields, report, logarithmic=False)
    _images(output, fields, report, logarithmic=True)
    hemisphere = report['far_hemisphere']
    exposure_s = report['exposure_ms']*1e-3
    lines = [
        'Independent component comparison',
        'Particle only: sphere above a smooth wafer.',
        'Haze only: saved rough wafer without a sphere.',
        'These are separate experiments. Do not add their powers to reconstruct the combined case.',
        '',
        'Upper hemisphere (all directions):',
        '  Particle DSC: {:.8g} um^2; power: {:.8g} pW; energy: {:.8g} fJ'.format(
            hemisphere['particle_dsc_um2'], hemisphere['particle_equivalent_power_w']*1e12,
            hemisphere['particle_equivalent_power_w']*exposure_s*1e15),
        '  Haze DSC:     {:.8g} um^2; power: {:.8g} pW; energy: {:.8g} fJ'.format(
            hemisphere['haze_dsc_um2'], hemisphere['haze_equivalent_power_w']*1e12,
            hemisphere['haze_equivalent_power_w']*exposure_s*1e15),
        '',
        'Full image-plane integrated power (independent objective channels):',
        '  Channel    NA     Tilt      Particle (pW)    Haze (pW)    Particle (fJ)    Haze (fJ)']
    for channel in report['channels']:
        lines.append('  Ch{channel}       {na:.4f}  {tilt_deg:+6.1f}       '
                     '{particle_power_w:.8g}       {haze_power_w:.8g}       '
                     '{particle_energy_j:.8g}       {haze_energy_j:.8g}'.format(
                         **dict(channel, particle_power_w=channel['particle_power_w']*1e12,
                                haze_power_w=channel['haze_power_w']*1e12,
                                particle_energy_j=channel['particle_energy_j']*1e15,
                                haze_energy_j=channel['haze_energy_j']*1e15)))
    lines += ['', 'For an exposure of {:.6g} ms, energy (J) = power (W) x exposure (s).'.format(
        report['exposure_ms']),
        'Channels may overlap in angle; do not sum their powers as a unique collected total.',
        'Particle and haze have separate numerical color scales so weak patterns remain visible.',
        'Image channels share a scale within each component; use colorbar numbers or totals for cross-component comparisons.',
        'Only scattered light is modeled; specular leakage and physical convergence are not established.']
    (output/'totals.txt').write_text('\n'.join(lines)+'\n', encoding='utf-8')
