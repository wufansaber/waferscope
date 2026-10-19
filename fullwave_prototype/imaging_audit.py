"""Audit saved FEM images without using production geometry, lens rotations or FFT."""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np

from .direct_imaging import direct_image, polar_cap, ray_frame, rotate_between


def relative(actual, reference):
    if not np.isfinite(actual).all() or not np.isfinite(reference).all():
        raise ValueError('Non-finite field in independent comparison')
    return float(np.linalg.norm(actual-reference)/max(np.linalg.norm(reference), 1e-300))


def check_hankel():
    from .theory_optics import hankel_slice
    uv, weights, _ = polar_cap(.8, 0, 64, 192)
    local, _, s, p, _ = ray_frame(uv, 0)
    z = np.broadcast_to((0., 0., 1.), local.shape)
    far = rotate_between(np.broadcast_to((1., 0., 0.), local.shape), z, local)
    far *= np.sqrt(local[:, 2])[:, None]
    amplitude = np.column_stack((np.sum(far*s, axis=1), np.sum(far*p, axis=1)))
    reports = []
    for m, dz in ((1, 0), (1, .4), (20, 0)):
        x = np.linspace(-1.5, 1.5, 31)*m
        e, h, _ = direct_image(amplitude, uv, weights, np.column_stack((x, np.zeros_like(x))),
                              .355, magnification=m, defocus_um=dz)
        expected_e, expected_h = hankel_slice(x, .8, .355, m, dz)
        reports.append(dict(magnification=m, defocus_um=dz,
                            e_relative_l2=relative(e, expected_e), h_relative_l2=relative(h, expected_h)))
    return reports


def check_saved(source):
    source = Path(source)
    report = json.loads((source/'imaging_report.json').read_text(encoding='utf-8'))
    settings = report['settings']
    p = report['scattering_parameters']
    results, arrays = [], {}
    with np.load(source/'image_fields.npz') as data:
        for channel, tilt in enumerate(settings['tilt_degrees'], 1):
            prefix = 'ch{}_'.format(channel)
            axis = data[prefix+'image_axis_um']
            indices = np.flatnonzero(abs(axis) <= 12)
            xx, yy = np.meshgrid(axis[indices], axis[indices])
            xy = np.column_stack((xx.ravel(), yy.ravel()))
            u = data[prefix+'pupil_u_axis']
            uu, vv = np.meshgrid(u, u)
            mask = data[prefix+'pupil_mask']
            uv = np.column_stack((uu[mask], vv[mask]))
            weights = np.full(len(uv), (u[1]-u[0])**2)
            # Independent analytic OPD, not the production wavefront helper.
            na = settings['numerical_apertures'][channel-1]
            a, c, spherical = settings.get('aberrations_nm', (0, 0, 0))
            x, y = (uv/max(na, 1e-30)).T
            r2 = x*x+y*y
            opd = a*np.sqrt(6)*(x*x-y*y)+c*np.sqrt(8)*x*(3*r2-2)+spherical*np.sqrt(5)*(6*r2*r2-6*r2+1)
            phase = np.exp(2j*np.pi*opd/p['wavelength_nm'])
            row = dict(channel=channel, tilt_deg=tilt, samples=len(xy))
            row['complex_fields_checked'] = all(prefix+'image_'+component+'_'+name in data
                for component in ('e', 'hz0') for name in ('with', 'without'))
            fields = {}
            for name in ('with', 'without'):
                e, h, flux = direct_image(data[prefix+'amplitude_'+name][mask]*phase[:, None], uv, weights,
                    xy, p['wavelength_nm']/1000, tilt, settings['magnification'],
                    settings['image_defocus_um'], settings['transmission'], settings['analyzer'])
                ix = np.ix_(indices, indices)
                fft_flux = data[prefix+'image_'+name+'_w_um2'][ix]
                direct = flux.reshape(xx.shape)*settings['irradiance_w_um2']
                row[name] = dict(flux_relative_l2=relative(direct, fft_flux))
                if row['complex_fields_checked']:
                    fft_e = data[prefix+'image_e_'+name][ix].reshape(-1, 3)
                    fft_h = data[prefix+'image_hz0_'+name][ix].reshape(-1, 3)
                    row[name].update(e_relative_l2=relative(e, fft_e), h_relative_l2=relative(h, fft_h))
                arrays[prefix+name+'_direct'] = direct
                arrays[prefix+name+'_fft'] = fft_flux
                fields[name] = e, h
            delta_e = fields['with'][0]-fields['without'][0]
            delta_h = fields['with'][1]-fields['without'][1]
            arrays[prefix+'change_direct'] = np.cross(delta_e, delta_h.conj())[:, 2].real.reshape(xx.shape)*settings['irradiance_w_um2']
            arrays[prefix+'axis_um'] = axis[indices]
            mid = len(indices)//2
            positive = np.flatnonzero((axis[indices] > 0) & (axis[indices] <= 6))
            total = arrays[prefix+'with_direct']
            if len(positive) and np.max(total) > 0:
                dark = positive[np.argmin(total[positive, mid])]
                row['positive_y_dip'] = dict(x_um=float(axis[indices][mid]),
                    y_um=float(axis[indices][dark]), flux_w_um2=float(total[dark, mid]),
                    fraction_of_roi_peak=float(total[dark, mid]/total.max()))
            row['passed'] = max(v for name in ('with', 'without') for v in row[name].values()) < 1e-10
            results.append(row)
    return dict(source=str(source.resolve()), channels=results,
        scope='Same sampled FEM pupil; independent coordinate rotations and direct quadrature. '
              'No fitted phase, amplitude normalization, FFT or image alignment. '
              'Compact archives check flux only, marked by complex_fields_checked=false. '
              'This tests imaging implementation, not convergence of the FEM source.'), arrays


def bv_reference(p, weights, directions, output, python, lmax=14):
    if p.get('defect_shape', 'sphere') != 'sphere':
        raise ValueError('The Bobbert-Vlieger image reference requires a sphere')
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output/'directions.npz', directions=directions, weights=weights)
    request = dict(parameters=p, kind='bv', lmax=lmax,
                   directions_file=str((output/'directions.npz').resolve()))
    (output/'request.json').write_text(json.dumps(request), encoding='utf-8')
    with (output/'reference.log').open('w', encoding='utf-8') as log:
        subprocess.run([python, '-m', 'fullwave_prototype.audit_reference', '--request',
                        str(output/'request.json'), '--output', str(output/'reference.npz')],
                       check=True, timeout=600, stdout=log, stderr=subprocess.STDOUT)
    with np.load(output/'reference.npz') as saved:
        amplitude = saved['amplitude']
    height = (p['diameter_nm']/2+p['gap_nm'])/1000
    return amplitude*np.exp(-2j*np.pi/(p['wavelength_nm']/1000)*height*directions[:, 2])[:, None]


def check_flat(source, output, python):
    source, output = Path(source), Path(output)
    results, arrays = [], {}
    for pol in ('p', 's'):
        folder = source/('flat_'+pol)
        saved = json.loads((folder/'report.json').read_text(encoding='utf-8'))
        p = saved['parameters']
        with np.load(folder/'images.npz') as data:
            axis = data['axis_um']
            indices = np.flatnonzero(abs(axis) <= 12)
            x, y = np.meshgrid(axis[indices], axis[indices])
            xy = np.column_stack((x.ravel(), y.ravel()))
            arrays['flat_axis_um'] = axis[indices]
            for channel, tilt in enumerate((45, -45, 0), 1):
                print('Independent BV image: {} polarization, channel {}'.format(pol, channel), flush=True)
                prefix = '{}_ch{}'.format(pol, channel)
                quadratures = []
                for order, azimuths, lmax in ((24, 96, 14), (48, 192, 14), (48, 192, 18)):
                    uv, weights, directions = polar_cap(.5, tilt, order, azimuths)
                    amplitude = bv_reference(p, weights, directions,
                        output/(prefix+'_q{}_l{}'.format(order, lmax)), python, lmax)
                    _, _, flux = direct_image(amplitude, uv, weights, xy, p['wavelength_nm']/1000, tilt)
                    quadratures.append(flux.reshape(x.shape))
                reference = quadratures[-1]
                fine_axis = np.linspace(-12, 12, 49)
                fx, fy = np.meshgrid(fine_axis, fine_axis)
                _, _, fine = direct_image(amplitude, uv, weights,
                    np.column_stack((fx.ravel(), fy.ravel())), p['wavelength_nm']/1000, tilt)
                fine = fine.reshape(fx.shape)
                fem = data['ch{}_fem'.format(channel)][np.ix_(indices, indices)]
                bv = data['ch{}_bv'.format(channel)][np.ix_(indices, indices)]
                row = dict(polarization=pol, channel=channel, tilt_deg=tilt,
                    fft_pupil_samples=saved.get('imaging_settings', {}).get('pupil_samples', 32),
                    quadrature_relative_l2=relative(quadratures[0], quadratures[1]),
                    multipole_relative_l2=relative(quadratures[1], reference),
                    fem_fft_vs_bv_direct_relative_l2=relative(fem, reference),
                    bv_fft_vs_bv_direct_relative_l2=relative(bv, reference))
                row['reference_converged'] = max(row['quadrature_relative_l2'], row['multipole_relative_l2']) < 1e-4
                row['fem_agrees_with_reference'] = row['fem_fft_vs_bv_direct_relative_l2'] < .02
                mid = len(fine_axis)//2
                positive = np.flatnonzero((fine_axis > 0) & (fine_axis <= 6))
                dark = positive[np.argmin(fine[positive, mid])]
                row['positive_y_dip'] = dict(x_um=0., y_um=float(fine_axis[dark]),
                    flux_w_um2=float(fine[dark, mid]*1e-6), fraction_of_roi_peak=float(fine[dark, mid]/fine.max()))
                results.append(row)
                arrays[prefix+'_direct'] = reference*1e-6
                arrays[prefix+'_fine_direct'] = fine*1e-6
                arrays['flat_fine_axis_um'] = fine_axis
                arrays[prefix+'_fem'] = fem*1e-6
    return dict(channels=results, irradiance_w_um2=1e-6,
        scope='Flat wafer only. Independent SCATMECH BV amplitudes and polar-Gauss direct imaging '
              'versus saved FEM/FFT; finite-pupil-grid error remains in the latter. '
              'Sphere-origin translation is known, not fitted. Absolute flux; no spatial registration.'), arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, help='Saved rough imaging directory')
    parser.add_argument('--flat-source', required=True, help='Directory containing flat_p and flat_s')
    parser.add_argument('--reference-python', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    hankel = check_hankel()
    saved, arrays = check_saved(args.source)
    flat, other = check_flat(args.flat_source, output/'bv', args.reference_python)
    arrays.update(other)
    report = dict(hankel=hankel, saved=saved, flat=flat)
    np.savez_compressed(output/'fields.npz', **arrays)
    (output/'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    from .imaging_audit_plot import plot
    plot(output, arrays)
    print(json.dumps(report, indent=2), flush=True)
    if (not all(r['passed'] for r in saved['channels'])
            or not all(r['reference_converged'] and r['fem_agrees_with_reference'] for r in flat['channels'])
            or max(max(r['e_relative_l2'], r['h_relative_l2']) for r in hankel) >= 1e-10):
        raise RuntimeError('Independent image audit failed; inspect report.json')


if __name__ == '__main__':
    main()
