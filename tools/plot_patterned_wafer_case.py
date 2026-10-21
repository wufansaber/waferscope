"""Absolute-scale sample/reference views for the finite etched-line example."""
import json

import numpy as np


def plot_case(output, data, arrays, report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, TwoSlopeNorm
    from matplotlib.patches import Rectangle, Circle
    pattern = report['pattern']
    physical, settings = report['scattering_parameters'], report['settings']
    diameter_nm = physical['diameter_nm']
    half = pattern['half_width_um']
    warning = ('NOT CONVERGED' if report.get('mesh_validation', {}).get('status') == 'not_converged'
               else 'UNVALIDATED FEM PREVIEW')
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.6), layout='constrained')
    for c in (-pattern['pitch_um']/2, pattern['pitch_um']/2):
        axes[0].add_patch(Rectangle((-pattern['line_length_um']/2, c-pattern['trench_width_um']/2),
            pattern['line_length_um'], pattern['trench_width_um'], facecolor='#9aabb8', edgecolor='black'))
    axes[0].add_patch(Circle((0, 0), diameter_nm/2000, color='#c43f75'))
    axes[0].set(xlim=(-half, half), ylim=(-half, half), aspect='equal', xlabel='Object x (um)',
                ylabel='Object y (um)', title='Assumed geometry | two finite trenches')
    axes[0].text(-.94*half, .92*half,
        'Pitch {:g} nm; line {:g} nm; depth {:g} nm\n{:g} nm sphere on central ridge'.format(
            pattern['pitch_um']*1000, (pattern['pitch_um']-pattern['trench_width_um'])*1000,
            pattern['trench_depth_um']*1000, diameter_nm),
        fontsize=8, va='top', bbox=dict(facecolor='white', alpha=.85, edgecolor='none'))
    limit = max(np.nanmax(data['disk_'+key]) for key in ('total', 'background'))
    for ax, key, title in zip(axes[1:], ('background', 'total'), ('Lines only', 'Lines + {:g} nm sphere'.format(diameter_nm))):
        artist = ax.pcolormesh(data['disk_axis'], data['disk_axis'], data['disk_'+key],
                              norm=Normalize(0, limit), cmap='viridis', shading='auto')
        ax.set(xlabel='kx/k', ylabel='ky/k', aspect='equal', title=title+' | far field')
    fig.colorbar(artist, ax=axes[1:].tolist(), label='DSC (um^2/sr), common linear scale', shrink=.8)
    fig.suptitle('{} | {} nm, {} deg {} | No random haze | Not SEM reconstruction'.format(
        warning, physical['wavelength_nm'], physical['incidence_deg'], physical['polarization']), color='#a12020')
    fig.savefig(output/'geometry_far_field.png', dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(3, 3, figsize=(12, 11), layout='constrained')
    common = max(arrays['ch{}_image_{}_w_um2'.format(i, state)].max()
                 for i in (1, 2, 3) for state in ('with', 'without'))
    difference = max(abs(arrays['ch{}_image_difference_w_um2'.format(i)]).max() for i in (1, 2, 3))
    for row, (key, title) in enumerate((('without', 'Lines only'), ('with', 'Lines + sphere'),
                                       ('difference', 'Sample - reference (signed)'))):
        for i, ax in enumerate(axes[row], 1):
            axis = arrays['ch{}_image_axis_um'.format(i)]
            view_limit = settings['detector']['pixels_x']*settings['detector']['pixel_pitch_um']/2
            selected = abs(axis) <= view_limit
            value = arrays['ch{}_image_{}_w_um2'.format(i, key)][np.ix_(selected, selected)]
            norm = TwoSlopeNorm(vmin=-difference, vcenter=0, vmax=difference) if row == 2 else Normalize(0, common)
            artist = ax.pcolormesh(axis[selected], axis[selected], value, shading='auto',
                norm=norm, cmap='RdBu_r' if row == 2 else 'magma')
            ax.set(xlabel='Image x (um)', ylabel='Image y (um)', aspect='equal',
                   title='{} | CH{}'.format(title, i))
        fig.colorbar(artist, ax=axes[row].tolist(), label='W/um^2, common absolute scale', shrink=.8)
    fig.suptitle('Pitch {:g} nm | line {:g} nm | particle {:g} nm | {:g}x | NA {} | no random haze\n'.format(
        pattern['pitch_um']*1000, (pattern['pitch_um']-pattern['trench_width_um'])*1000, diameter_nm,
        settings['magnification'], '/'.join('{:g}'.format(na) for na in settings['numerical_apertures']))+
        warning+': absolute brightness and detectability are not validated', color='#a12020')
    fig.savefig(output/'optical_comparison.png', dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(12, 7.5), layout='constrained')
    low = min(arrays['ch{}_detector_gray_first_{}_dn'.format(i, state)].min()
              for i in (1, 2, 3) for state in ('sample', 'reference'))
    high = max(arrays['ch{}_detector_gray_first_{}_dn'.format(i, state)].max()
               for i in (1, 2, 3) for state in ('sample', 'reference'))
    for row, state in enumerate(('reference', 'sample')):
        for i, ax in enumerate(axes[row], 1):
            value = arrays['ch{}_detector_gray_first_{}_dn'.format(i, state)]
            detector = settings['detector']
            lx, ly = (detector['pixels_x']*detector['pixel_pitch_um']/2,
                      detector['pixels_y']*detector['pixel_pitch_um']/2)
            artist = ax.imshow(value, origin='lower', extent=(-lx, lx, -ly, ly),
                cmap='gray', vmin=low, vmax=high, interpolation='nearest')
            ax.set(title='{} raw DN | CH{}'.format(state.capitalize(), i),
                   xlabel='Sensor x (um)', ylabel='Sensor y (um)')
    fig.colorbar(artist, ax=axes.ravel().tolist(), label='Raw DN (same display window)', shrink=.85)
    exposure_us = report['settings']['detector']['exposure_ms']*1000
    fig.suptitle('Gemi area approximation | {:g} um pixels, {}x{} ROI, {:g} us | {:g}x | no TDI\n'.format(
        detector['pixel_pitch_um'], detector['pixels_x'], detector['pixels_y'], exposure_us, settings['magnification'])+
        warning+': gray frames and SNR are not validated physical predictions', color='#a12020')
    fig.savefig(output/'camera_comparison.png', dpi=180)
    plt.close(fig)
    rows = []
    for i, ch in enumerate(report['channels'], 1):
        rows.append(dict(channel=i, background_power_w=ch['without']['image_power_w'],
            sample_power_w=ch['with']['image_power_w'],
            signed_power_difference_w=ch['with']['image_power_w']-ch['without']['image_power_w'],
            gray_snr=ch['detector']['gray_snr'], gray_valid=ch['detector']['gray_valid'],
            decomposition_relative_error=ch['decomposition_relative_max_error']))
    (output/'case_summary.json').write_text(json.dumps(dict(status=warning, parameters=physical,
        pattern=pattern, settings=settings, channels=rows,
        warning='Synthetic finite pattern; not reconstructed SEM or a convergence-certified result.'), indent=2), encoding='utf-8')
