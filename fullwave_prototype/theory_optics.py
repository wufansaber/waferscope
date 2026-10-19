"""PEC image-dipole and nonparaxial circular-pupil Hankel benchmarks."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import jv

from optical_propagation import pupil_geometry, propagate_pupil
from .far_field import radiation_amplitude
from .imaging import ImagingSettings
from .parameters import Parameters


def pec_image_check():
    theta = np.linspace(.02, 1.4, 50)
    phi = np.linspace(0, 2*np.pi, 36, endpoint=False)
    t, f = np.meshgrid(theta, phi)
    k = np.stack((np.sin(t)*np.cos(f), np.sin(t)*np.sin(f), np.cos(t)), -1).reshape(-1, 3)
    s = np.column_stack((-np.sin(f.ravel()), np.cos(f.ravel()), np.zeros(f.size)))
    basis_p = np.cross(k, s)
    position = np.array((.015, -.02, .04))
    mirror = position*np.array((1, 1, -1))
    results = []
    for moment in (np.array((1, 0, 0)), np.array((0, 0, 1))):
        image_moment = moment*np.array((-1, -1, 1))
        for loss in (1e3, 1e5):
            p = Parameters(substrate_n=1, substrate_k=loss)
            actual = radiation_amplitude(position[None], moment[None], k, p)
            source = (np.exp(-1j*p.k0*(k@position))[:, None]*moment+
                      np.exp(-1j*p.k0*(k@mirror))[:, None]*image_moment)
            exact = p.k0**2/(4*np.pi)*np.column_stack((np.sum(s*source, 1), np.sum(basis_p*source, 1)))
            results.append(dict(orientation='horizontal' if moment[0] else 'vertical',
                substrate_k=loss, relative_complex_l2=float(np.linalg.norm(actual-exact)/np.linalg.norm(exact))))
    return results


def hankel_slice(x, na, wavelength, magnification, defocus, order=200):
    """Analytic azimuthal integration of a uniform x-polarized collimated disk."""
    nodes, weights = leggauss(order)
    r, dr = (nodes+1)*na/2, weights*na/2
    s = r/magnification
    w = np.sqrt(1-s*s)
    argument = 2*np.pi/wavelength*np.outer(x, s)
    j0, j1, j2 = (jv(i, argument) for i in (0, 1, 2))
    weight = 2*np.pi/(wavelength*magnification)*r*dr/np.sqrt(w)
    weight = weight*np.exp(2j*np.pi/wavelength*defocus*(w-1))
    ex = (j0-s*s/(2*(1+w))*(j0-j2))@weight
    ez = (-1j*s*j1)@weight
    hy = (j0-s*s/(2*(1+w))*(j0+j2))@weight
    return np.column_stack((ex, np.zeros_like(ex), ez)), np.column_stack((np.zeros_like(ex), hy, np.zeros_like(ex)))


def vector_focus_check(samples=128):
    reports, arrays = [], {}
    for m, defocus in ((1., 0.), (1., .4), (20., 0.)):
        settings = ImagingSettings(magnification=m, pupil_samples=samples, padding=8, image_defocus_um=defocus)
        local, _, _, mask, du = pupil_geometry(.8, 0, samples)
        pupil = np.zeros(mask.shape+(2, 2), complex)
        pupil[..., 0, 0] = mask
        axis, e, h = propagate_pupil(pupil, local, du, .355, settings)
        roi = abs(axis) <= 2*m
        x = axis[roi]
        actual_e, actual_h = e[len(axis)//2, roi, :, 0], h[len(axis)//2, roi, :, 0]
        expected_e, expected_h = hankel_slice(x, .8, .355, m, defocus)
        checked_e, checked_h = hankel_slice(x, .8, .355, m, defocus, 320)
        flux = (actual_e[:, 0]*actual_h[:, 1].conj()).real
        exact_flux = (expected_e[:, 0]*expected_h[:, 1].conj()).real
        reports.append(dict(magnification=m, image_na=.8/m, defocus_um=defocus,
            electric_relative_l2=float(np.linalg.norm(actual_e-expected_e)/np.linalg.norm(expected_e)),
            magnetic_relative_l2=float(np.linalg.norm(actual_h-expected_h)/np.linalg.norm(expected_h)),
            flux_relative_l2=float(np.linalg.norm(flux-exact_flux)/np.linalg.norm(exact_flux)),
            oracle_quadrature_relative=float(np.linalg.norm(checked_e-expected_e)/np.linalg.norm(checked_e))))
        key = 'm{}_z{}'.format(m, defocus)
        arrays[key+'_axis'], arrays[key+'_flux'], arrays[key+'_reference'] = x, flux, exact_flux
    return reports, arrays
