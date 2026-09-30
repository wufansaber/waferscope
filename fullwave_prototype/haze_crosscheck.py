"""Compare saved FEM haze with NIST roughness models on one angular grid."""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from .angular_plotting import angular_disk_edges, format_angular_disk


def _surface_transform(triangles, theta, phi, wavelength_um, incidence_deg, order):
    a, b, c = (triangles[:, index, :] for index in range(3))
    determinant = np.abs((b[:, 0]-a[:, 0])*(c[:, 1]-a[:, 1]) -
                         (b[:, 1]-a[:, 1])*(c[:, 0]-a[:, 0]))
    nodes, weights = leggauss(order)
    nodes, weights = (nodes+1)/2, weights/2
    positions, weighted_heights = [], []
    area = mean = second = 0.
    for u, wu in zip(nodes, weights):
        for v, wv in zip(nodes, weights):
            point = a+u*(b-a)+(1-u)*v*(c-a)
            weight = determinant*(1-u)*wu*wv
            positions.append(point[:, :2])
            weighted_heights.append(weight*point[:, 2])
            area += weight.sum()
            mean += np.sum(weight*point[:, 2])
            second += np.sum(weight*point[:, 2]**2)
    positions = np.concatenate(positions)
    weighted_heights = np.concatenate(weighted_heights)
    polar, azimuth = np.meshgrid(theta, phi, indexing='ij')
    frequency = np.column_stack((np.sin(polar).ravel()*np.cos(azimuth).ravel(),
        np.sin(polar).ravel()*np.sin(azimuth).ravel()-np.sin(np.deg2rad(incidence_deg))))/wavelength_um
    transform = np.empty(len(frequency), complex)
    for start in range(0, len(frequency), 128):
        transform[start:start+128] = (np.exp(-2j*np.pi*(frequency[start:start+128] @ positions.T))
                                      @ weighted_heights)
    return transform.reshape(polar.shape), dict(area_um2=float(area),
        mean_nm=float(mean/area*1000), rms_nm=float(np.sqrt(second/area)*1000))


def _nist_dsc(theta, phi, parameters, spectrum, area_um2):
    from pySCATMECH.brdf import BRDF_Model
    from pySCATMECH.mueller import Polarization

    model = BRDF_Model('Microroughness_BRDF_Model', {
        'lambda': parameters['wavelength_nm']/1000,
        'substrate': complex(parameters['substrate_n'], parameters['substrate_k']),
        'type': 0, 'psd': spectrum})
    incident = np.deg2rad(parameters['incidence_deg'])
    polarization = Polarization(parameters['polarization'])
    result = np.zeros((len(theta), len(phi)))
    for i, polar in enumerate(theta):
        if np.cos(polar) <= 1e-12:
            continue
        for j, azimuth in enumerate(phi):
            # SCATMECH's incidence-plane X is global +Y for illumination from -KY.
            local_azimuth = (azimuth-np.pi/2) % (2*np.pi)
            brdf = float((model.MuellerBRDF(incident, polar, local_azimuth,
                                           coords='psps') @ polarization)[0])
            result[i, j] = max(brdf, 0.)*area_um2*np.cos(incident)*np.cos(polar)
    return result


def _hemisphere(values, theta):
    return float(2*np.pi*np.trapz(np.mean(values, axis=1)*np.sin(theta), theta))


def _channel_integrals(values, theta, phi, apertures, tilts):
    polar, azimuth = np.meshgrid(theta, phi, indexing='ij')
    x = np.sin(polar)*np.cos(azimuth)
    z = np.cos(polar)
    result = []
    for na, tilt in zip(apertures, tilts):
        if na == 0:
            result.append(0.)
            continue
        angle = np.deg2rad(tilt)
        mask = x*np.sin(angle)+z*np.cos(angle) >= np.sqrt(1-na*na)
        result.append(_hemisphere(values*mask, theta))
    return result


def _peak(values, theta, phi):
    index = np.unravel_index(np.argmax(values), values.shape)
    return dict(theta_deg=float(np.rad2deg(theta[index[0]])),
                phi_deg=float(np.rad2deg(phi[index[1]])),
                dsc_um2_sr=float(values[index]))


def _similarity(first, second, theta):
    product = _hemisphere(first*second, theta)
    return float(product/np.sqrt(_hemisphere(first*first, theta)*
                                 _hemisphere(second*second, theta)))


def _plot(output, theta, phi, maps, metrics, parameters):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, PowerNorm

    names = ('FEM: one tapered surface', 'NIST: same height map',
             'NIST: Gaussian ensemble PSD', 'NIST: ABC PSD')
    keys = ('fem', 'same_height', 'gaussian', 'abc')
    maximum = max(float(np.max(values)) for values in maps.values())
    x_edges, y_edges = angular_disk_edges(theta, phi)
    for shared_log in (False, True):
        fig, axes = plt.subplots(2, 2, figsize=(17 if not shared_log else 13, 11),
                                 layout='constrained')
        for ax, key, name in zip(axes.ravel(), keys, names):
            local_max = max(float(np.max(maps[key])), 1e-30)
            norm = (LogNorm(vmin=max(maximum*1e-6, 1e-30), vmax=max(maximum, 1e-30))
                    if shared_log else PowerNorm(.5, vmin=0, vmax=local_max))
            image = ax.pcolormesh(x_edges, y_edges, maps[key],
                                  norm=norm, cmap='viridis', shading='flat', rasterized=True)
            format_angular_disk(ax)
            ax.scatter([0], [-np.sin(np.deg2rad(parameters['incidence_deg']))],
                       s=65, marker='*', facecolor='white', edgecolor='black', zorder=8)
            ax.set_title('{}\nHemisphere DSC = {:.4g} um^2'.format(
                name, metrics[key]['hemisphere_dsc_um2']))
            if not shared_log:
                fig.colorbar(image, ax=ax, shrink=.65, label='DSC (um^2 / sr)')
        if shared_log:
            fig.colorbar(image, ax=axes.ravel().tolist(), shrink=.72,
                         label='DSC (um^2 / sr), one shared logarithmic scale')
        fig.suptitle('Haze only | global kx/k, ky/k | source -KY | {:g} nm | {:g} deg {} | {}'.format(
            parameters['wavelength_nm'], parameters['incidence_deg'], parameters['polarization'],
            'shared logarithmic scale' if shared_log else 'separate numerical color scales'), fontsize=13)
        fig.savefig(output/('haze_angular_shared_log.png' if shared_log else
                            'haze_angular_comparison.png'), dpi=150)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5), layout='constrained')
    x = np.arange(3)
    width = .19
    for index, (key, name) in enumerate(zip(keys, names)):
        ax.bar(x+(index-1.5)*width, metrics[key]['channels_dsc_um2'], width,
               label=name)
    ax.set(xticks=x, xticklabels=('Ch1 +45 deg', 'Ch2 -45 deg', 'Ch3 0 deg'),
           ylabel='Collected haze DSC (um^2), 5-deg angular-grid estimate')
    ax.legend(fontsize=8)
    ax.grid(axis='y', alpha=.2)
    fig.savefig(output/'haze_channel_comparison.png', dpi=150)
    plt.close(fig)


def compare(source, output, spot_length_um=2., spot_width_um=2.):
    source, output = Path(source), Path(output)
    report = json.loads((source/'report.json').read_text(encoding='utf-8'))
    parameters, roughness = report['parameters'], report['roughness']
    if parameters['beam_waist_um'] != 0 or roughness['model'] not in ('legacy', 'gaussian'):
        raise ValueError('This cross-check requires plane-wave illumination and a saved rough surface')
    with np.load(source/'paired_fields.npz', allow_pickle=False) as saved:
        theta, phi = saved['theta'], saved['phi']
        maps = {'fem': saved['dsc_background']}
        triangles = saved['surface_triangles_um']
    wavelength_um = parameters['wavelength_nm']/1000
    transform, surface = _surface_transform(triangles, theta, phi, wavelength_um,
        parameters['incidence_deg'], order=6)
    lower_order, _ = _surface_transform(triangles, theta, phi, wavelength_um,
        parameters['incidence_deg'], order=4)
    area = surface['area_um2']
    if not np.isclose(area, roughness['width_um']**2, rtol=1e-10):
        raise ValueError('FEM height-map area differs from reported rough patch width')
    unit = _nist_dsc(theta, phi, parameters, {None: 'Unit_PSD_Function'}, area)
    maps['same_height'] = unit*(abs(transform)**2/area)
    sigma = roughness['rms_nm']/1000
    length = roughness['correlation_length_um']
    maps['gaussian'] = _nist_dsc(theta, phi, parameters, {
        None: 'Gaussian_PSD_Function', 'sigma': sigma, 'length': length}, area)
    maps['abc'] = _nist_dsc(theta, phi, parameters, {
        None: 'ABC_PSD_Function', 'A': sigma**2*length**2/(2*np.pi),
        'B': length, 'C': 3.}, area)
    image_report = json.loads((source.parent/'imaging_report.json').read_text(encoding='utf-8'))
    settings = image_report['settings']
    metrics = {key: dict(hemisphere_dsc_um2=_hemisphere(values, theta),
                         peak=_peak(values, theta, phi),
                         channels_dsc_um2=_channel_integrals(values, theta, phi,
                             settings['numerical_apertures'], settings['tilt_degrees']))
               for key, values in maps.items()}
    metrics['same_height']['map_similarity_to_fem'] = _similarity(
        maps['fem'], maps['same_height'], theta)
    metrics['gaussian']['map_similarity_to_fem'] = _similarity(
        maps['fem'], maps['gaussian'], theta)
    metrics['abc']['map_similarity_to_fem'] = _similarity(
        maps['fem'], maps['abc'], theta)
    result = dict(source=str(source.resolve()), parameters={key: parameters[key] for key in
        ('wavelength_nm', 'incidence_deg', 'polarization', 'substrate_n', 'substrate_k')},
        roughness={key: roughness[key] for key in ('model', 'rms_nm', 'correlation_length_um',
                                                   'width_um', 'grid_points', 'seed')},
        surface=surface, gaussian_spot_area_um2=spot_length_um*spot_width_um,
        surface_transform_relative_l2=float(np.linalg.norm(transform-lower_order)/
                                            np.linalg.norm(transform)),
        map_similarity_definition='Solid-angle-weighted cosine similarity of nonnegative DSC maps; 1 means the same angular shape up to scale.',
        map_scope='Haze only. NIST maps use A=FEM geometric patch area. Unit-PSD map uses the exact FEM height realization. Statistical Gaussian L is not a legacy-FEM input.' if roughness['model'] == 'legacy' else 'Haze only. NIST maps use A=FEM geometric patch area. Unit-PSD map uses the exact FEM height realization.',
        spot_area_note='Current UI spot controls power scaling, not FEM rough-surface extent.',
        nist_gaussian_spot_area_hemisphere_dsc_um2=metrics['gaussian']['hemisphere_dsc_um2']*
            spot_length_um*spot_width_um/area,
        nist_abc_spot_area_hemisphere_dsc_um2=metrics['abc']['hemisphere_dsc_um2']*
            spot_length_um*spot_width_um/area,
        channels=dict(na=settings['numerical_apertures'], tilt_deg=settings['tilt_degrees'],
                      integration='5-degree angular-grid cone mask; not direct complex-pupil integration'),
        metrics=metrics)
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output/'haze_maps.npz', theta=theta, phi=phi, **maps)
    _plot(output, theta, phi, maps, metrics, parameters)
    (output/'haze_crosscheck.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    lines = ['# Haze-only FEM / NIST cross-check', '',
        'Both models use the saved FEM wavelength, Si index, incidence, polarization and 5-degree angular grid.',
        'NIST here is Microroughness_BRDF_Model only; no Bobbert-Vlieger particle is included.',
        'All four maps use the FEM rough patch geometric area ({:.6g} um^2).'.format(area),
        'FEM uses the {} roughness generator and one tapered height realization; Gaussian and ABC are ensemble PSD models.'.format(roughness['model']),
        'For legacy FEM, L is not used to generate the surface, so the statistical Gaussian and ABC columns are not matched physical inputs.' if roughness['model'] == 'legacy' else 'Gaussian FEM uses L for its stationary parent before finite-patch taper.',
        'The same-height NIST map uses the piecewise-planar FEM surface Fourier transform and Unit_PSD_Function.',
        '', '| Model | Hemisphere DSC (um^2) | Ratio to FEM | Peak (theta, phi) | Ch1 | Ch2 | Ch3 |',
        '| --- | ---: | ---: | --- | ---: | ---: | ---: |']
    for key, name in zip(('fem', 'same_height', 'gaussian', 'abc'),
                         ('FEM', 'NIST same height', 'NIST Gaussian PSD', 'NIST ABC PSD')):
        item = metrics[key]
        lines.append('| {} | {:.8g} | {:.4f} | ({:g}, {:g}) deg | {:.6g} | {:.6g} | {:.6g} |'.format(
            name, item['hemisphere_dsc_um2'],
            item['hemisphere_dsc_um2']/metrics['fem']['hemisphere_dsc_um2'],
            item['peak']['theta_deg'], item['peak']['phi_deg'], *item['channels_dsc_um2']))
    lines += ['', 'Saved FEM input RMS: {:.6g} nm; actual tapered patch RMS: {:.6g} nm.'.format(
        roughness['rms_nm'], surface['rms_nm']),
        'Surface Fourier quadrature relative L2 change (order 4 vs 6): {:.3g}.'.format(
            result['surface_transform_relative_l2']),
        'NIST Gaussian with {:.6g} x {:.6g} um spot area would yield hemisphere DSC {:.8g} um^2, '
        'but FEM does not simulate that full rough area.'.format(spot_length_um, spot_width_um,
            result['nist_gaussian_spot_area_hemisphere_dsc_um2']),
        'NIST ABC with that spot area would yield hemisphere DSC {:.8g} um^2.'.format(
            result['nist_abc_spot_area_hemisphere_dsc_um2']),
        'Angular shape similarity is the solid-angle-weighted cosine of two DSC maps (1 = identical shape up to scale).',
        'This is not an FEM mesh/PML convergence study or a statistical noise-variance validation.', '']
    (output/'README.md').write_text('\n'.join(lines), encoding='utf-8')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, help='Saved FEM scattering directory')
    parser.add_argument('--output', required=True)
    parser.add_argument('--spot-length-um', type=float, default=2.)
    parser.add_argument('--spot-width-um', type=float, default=2.)
    args = parser.parse_args()
    if args.spot_length_um <= 0 or args.spot_width_um <= 0:
        parser.error('Spot dimensions must be positive')
    result = compare(args.source, args.output, args.spot_length_um, args.spot_width_um)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
