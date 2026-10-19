"""Absolute image/collection checks against continuous diffraction quadratures."""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from optical_propagation import collimated_pupil, propagate_pupil, pupil_geometry
from .direct_imaging import direct_image, polar_cap, ray_frame
from .imaging import ImagingSettings, axial_flux
from .imaging_audit import check_saved, relative
from .theory_optics import hankel_slice


def dipole_amplitude(directions, moment=(.3, .2j, 1)):
    """Unit-normalized free-space electric dipole, without a substrate or fit."""
    directions = np.asarray(directions)
    moment = np.asarray(moment, complex)
    phi = np.arctan2(directions[:, 1], directions[:, 0])
    s = np.column_stack((-np.sin(phi), np.cos(phi), np.zeros(len(phi))))
    s[np.hypot(directions[:, 0], directions[:, 1]) < 1e-12] = (-1, 0, 0)
    p = np.cross(directions, s)
    return np.column_stack((s@moment, p@moment))


def axial_dipole_cap(na):
    """Exact integral of sin(theta)^2 over a normal circular spherical cap."""
    cosine = np.sqrt(1-na*na)
    height = na*na/(1+cosine)
    return 2*np.pi*height**2*(cosine+2)/3


def continuous_cap(na, tilt, order, azimuths):
    """Resolve the horizon analytically instead of masking discontinuous samples."""
    angle = np.deg2rad(tilt)
    boundary = np.cos(angle)
    if na <= boundary:
        return polar_cap(na, tilt, order, azimuths)
    nodes, weights = leggauss(order)
    radius = (nodes+1)*boundary/2
    dr = weights*boundary/2
    phi = np.arange(azimuths)*2*np.pi/azimuths
    uv_inner = np.stack((radius[:, None]*np.cos(phi), radius[:, None]*np.sin(phi)), -1).reshape(-1, 2)
    inner_weights = np.broadcast_to((radius*dr*2*np.pi/azimuths)[:, None],
                                    (order, azimuths)).ravel()
    # sin(t)^2 regularizes the square-root arc opening at the radial boundary.
    t = (nodes+1)*np.pi/4
    radius = boundary+(na-boundary)*np.sin(t)**2
    dr = weights*np.pi/4*(na-boundary)*np.sin(2*t)
    arc_start = np.arccos(boundary*np.sqrt(1-radius**2)/(abs(np.sin(angle))*radius))
    pn, pw = leggauss(azimuths)
    half_arc = np.pi-arc_start
    phi = np.pi+half_arc[:, None]*pn
    if tilt < 0:
        phi += np.pi
    uv_outer = np.stack((radius[:, None]*np.cos(phi), radius[:, None]*np.sin(phi)), -1).reshape(-1, 2)
    outer_weights = (radius*dr*half_arc)[:, None]*pw
    uv = np.concatenate((uv_inner, uv_outer))
    quadrature = np.concatenate((inner_weights, outer_weights.ravel()))
    _, directions, _, _, _ = ray_frame(uv, tilt)
    return uv, quadrature, directions


def dipole_case(na, tilt, samples, magnification=20, wavelength=.266):
    settings = ImagingSettings(magnification=magnification, pupil_samples=samples, padding=4)
    local, directions, rotation, mask, du = pupil_geometry(na, tilt, samples)
    jones = np.zeros(mask.shape+(2, 2), complex)
    jones[mask, :, 0] = dipole_amplitude(directions[mask])
    pupil = collimated_pupil(jones, directions, local, rotation)
    axis, e, h = propagate_pupil(pupil, local, du, wavelength, settings)
    flux = axial_flux(e[..., 0], h[..., 0])
    middle = len(axis)//2
    indices = np.flatnonzero(abs(axis) <= 1.5*wavelength*magnification/na + 1e-12)
    x = axis[indices]
    xy = np.concatenate((np.column_stack((x, np.zeros_like(x))),
                         np.column_stack((np.zeros_like(x), x))))
    actual = np.concatenate((flux[middle, indices], flux[indices, middle]))
    references, powers = [], []
    for order, azimuths in ((64, 128), (96, 192)):
        uv, weights, k = continuous_cap(na, tilt, order, azimuths)
        amplitude = dipole_amplitude(k)
        _, _, expected = direct_image(amplitude, uv, weights, xy, wavelength,
                                      tilt, magnification, chunk=16)
        references.append(expected)
        w = np.sqrt(1-np.sum(uv*uv, axis=1))
        powers.append(float(np.sum(np.sum(abs(amplitude)**2, axis=1)*weights/w)))
    collected = float(np.sum(abs(pupil)**2)*du**2)
    integrated = float(flux.sum()*(axis[1]-axis[0])**2)
    axial_grid = float(np.sum((1-directions[mask, 2]**2)/local[mask, 2])*du**2)
    row = dict(na=na, tilt_deg=tilt, samples=samples, magnification=magnification,
               image_na=na/magnification, wavelength_um=wavelength,
               collection_relative_error=abs(collected/powers[-1]-1),
               profile_relative_l2=relative(actual, references[-1]),
               reference_profile_relative_l2=relative(references[0], references[1]),
               reference_collection_relative_error=abs(powers[0]/powers[1]-1),
               parseval_relative_error=abs(integrated/collected-1))
    if tilt == 0:
        row['axial_dipole_exact_collection_relative_error'] = abs(axial_grid/axial_dipole_cap(na)-1)
    return row, dict(axis_um=x, fft=actual, reference=references[-1])


def uniform_disk_case(na, magnification, defocus, samples):
    settings = ImagingSettings(magnification=magnification, image_defocus_um=defocus,
                               pupil_samples=samples, padding=4)
    local, _, _, mask, du = pupil_geometry(na, 0, samples)
    pupil = np.zeros(mask.shape+(2, 2), complex)
    pupil[..., 0, 0] = mask
    axis, e, h = propagate_pupil(pupil, local, du, .266, settings)
    selected = abs(axis) <= 1.5*.266*magnification/na + 1e-12
    e, h = e[len(axis)//2, selected, :, 0], h[len(axis)//2, selected, :, 0]
    ref_e, ref_h = hankel_slice(axis[selected], na, .266, magnification, defocus, 240)
    fine_e, fine_h = hankel_slice(axis[selected], na, .266, magnification, defocus, 400)
    return dict(na=na, magnification=magnification, defocus_um=defocus, samples=samples,
                e_relative_l2=relative(e, ref_e), h_relative_l2=relative(h, ref_h),
                flux_relative_l2=relative(axial_flux(e, h), axial_flux(ref_e, ref_h)),
                reference_relative_l2=max(relative(ref_e, fine_e), relative(ref_h, fine_h)))


def plot(output, rows, arrays):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), constrained_layout=True)
    for na, tilt in ((.5, 0), (.96, 0), (.96, 45)):
        values = [r for r in rows if r['na'] == na and r['tilt_deg'] == tilt
                  and r['magnification'] == 20]
        label = 'NA={} / tilt={} deg'.format(na, tilt)
        for ax, key in zip(axes[:2], ('collection_relative_error', 'profile_relative_l2')):
            ax.plot([r['samples'] for r in values], [100*r[key] for r in values], 'o-', label=label)
    axes[0].set_ylabel('Collection error (%)')
    axes[1].set_ylabel('Image profile relative L2 (%)')
    for ax in axes[:2]:
        ax.set_xlabel('Pupil samples per axis')
        ax.set_xticks((32, 64, 96, 128))
        ax.grid(alpha=.25)
        ax.legend(fontsize=8)
    key = 'na0.96_tilt45_n128_m20'
    data = arrays[key]
    count = len(data['axis_um'])
    axes[2].plot(data['axis_um'], data['reference'][count:], label='Continuous direct integral')
    axes[2].plot(data['axis_um'], data['fft'][count:], 'o', ms=3, label='FFT, 128 samples')
    axes[2].set_xlabel('Image y (um), x=0')
    axes[2].set_ylabel('Flux / reference irradiance')
    axes[2].set_title('NA=0.96, tilt=45 deg, M=20')
    axes[2].legend(fontsize=8)
    axes[2].grid(alpha=.25)
    fig.savefig(output/'convergence.png', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--saved', action='append', default=[])
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    rows, arrays = [], {}
    for na, tilt, magnification in ((.5, 0, 20), (.5, 45, 20), (.5, -45, 20),
                                    (.96, 0, 20), (.96, 45, 20), (.96, 0, 1)):
        for samples in (32, 64, 96, 128):
            row, fields = dipole_case(na, tilt, samples, magnification)
            rows.append(row)
            key = 'na{}_tilt{}_n{}_m{}'.format(na, tilt, samples, magnification)
            arrays[key] = fields
            print(json.dumps(row), flush=True)
    disks = [uniform_disk_case(na, m, dz, n) for na, m, dz in
             ((.5, 20, 0), (.96, 20, 0), (.96, 1, 0), (.96, 1, .4))
             for n in (32, 64, 96, 128)]
    saved = [check_saved(source)[0] for source in args.saved]
    payload = dict(dipoles=rows, uniform_disks=disks, saved=saved,
                   scope='Absolute flux, no fitted normalization or phase. Profile checks use '
                         'two central line cuts, not a full 2D error bound. Tilted clipped caps '
                         'retain independently reported reference quadrature uncertainty.')
    (output/'report.json').write_text(json.dumps(payload, indent=2), encoding='utf-8')
    np.savez_compressed(output/'profiles.npz', **{key+'_'+name: value for key, fields
                       in arrays.items() for name, value in fields.items()})
    plot(output, rows, arrays)


if __name__ == '__main__':
    main()
