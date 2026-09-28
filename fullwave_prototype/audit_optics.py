"""Independent scalar Airy limit and near-field mesh comparison."""
import json
from pathlib import Path

import numpy as np
from scipy.special import j1

from optical_propagation import pupil_geometry, propagate_pupil
from .imaging import ImagingSettings, axial_flux


def airy_check():
    settings = ImagingSettings(pupil_samples=128, padding=8, magnification=100)
    na, wavelength = .2, .355
    local, _, _, mask, du = pupil_geometry(na, 0, settings.pupil_samples)
    pupil = np.zeros(mask.shape+(2, 2), complex)
    pupil[..., 0, 0] = mask
    axis, e, h = propagate_pupil(pupil, local, du, wavelength, settings)
    flux = axial_flux(e[..., 0], h[..., 0])
    row = flux[len(axis)//2]
    argument = 2*np.pi*na*axis/(wavelength*settings.magnification)
    pattern = np.ones_like(argument)
    nonzero = argument != 0
    pattern[nonzero] = (2*j1(argument[nonzero])/argument[nonzero])**2
    exact = (np.pi*na**2/(wavelength*settings.magnification))**2*pattern
    roi = abs(argument) < 12
    report = dict(na=na, magnification=settings.magnification, wavelength_um=wavelength,
                  absolute_slice_relative_l2=float(np.linalg.norm(row[roi]-exact[roi])/np.linalg.norm(exact[roi])),
                  peak_relative_error=float(abs(row.max()/exact.max()-1)),
                  note='Uniform circular collimated pupil, image NA=0.002; scalar Airy is the paraxial limit. Absolute amplitude retained, no fitted normalization.')
    return report, dict(axis_um=axis[roi], calculated=row[roi], airy=exact[roi])


def near_mesh_check(base, refined):
    with np.load(Path(base)/'paired_fields.npz') as f:
        a = dict(f)
    with np.load(Path(refined)/'paired_fields.npz') as f:
        b = dict(f)
    for key in ('near_y_um', 'near_z_um', 'surface_vertices_um'):
        np.testing.assert_allclose(a[key], b[key], rtol=0, atol=1e-14)
    pa = json.loads((Path(base)/'report.json').read_text())['parameters']
    pb = json.loads((Path(refined)/'report.json').read_text())['parameters']
    for key in ('diameter_nm', 'gap_nm', 'wavelength_nm', 'particle_n', 'particle_k',
                'substrate_n', 'substrate_k', 'incidence_deg', 'polarization'):
        if pa[key] != pb[key]:
            raise ValueError('Near-field cases differ in physics: '+key)
    yy, zz = np.meshgrid(a['near_y_um'], a['near_z_um'])
    radius = pa['diameter_nm']/2000
    center = radius+pa['gap_nm']/1000
    distance = np.sqrt(yy**2+(zz-center)**2)
    mask = (abs(zz) > .01) & (abs(distance-radius) > .005)
    x, y = a['near_e_with'], b['near_e_with']
    ia, ib = np.sum(abs(x)**2, axis=-1), np.sum(abs(y)**2, axis=-1)
    return dict(excluded_wafer_band_nm=10, excluded_sphere_band_nm=5, retained_points=int(mask.sum()),
                bulk_complex_e_relative_l2=float(np.linalg.norm(x[mask]-y[mask])/np.linalg.norm(y[mask])),
                bulk_intensity_relative_l2=float(np.linalg.norm(ia[mask]-ib[mask])/np.linalg.norm(ib[mask])),
                sampled_peak_base=float(ia.max()), sampled_peak_refined=float(ib.max()),
                sampled_peak_relative_change=float(abs(ia.max()/ib.max()-1)),
                bulk_sampled_peak_base=float(ia[mask].max()),
                bulk_sampled_peak_refined=float(ib[mask].max()),
                bulk_sampled_peak_relative_change=float(abs(ia[mask].max()/ib[mask].max()-1)),
                note='Legacy bulk mask uses distance to reference plane, not actual rough triangles. L2 and bulk_sampled_peak exclude these bands; sampled_peak is the UNMASKED YZ slice maximum, not a 3-D maximum. No interface/contact singularity validation.')


def main():
    output = Path('build/physical_audit')
    output.mkdir(parents=True, exist_ok=True)
    airy, arrays = airy_check()
    near = near_mesh_check('build/fem_roughness/order2', 'build/fem_roughness/refined')
    np.savez_compressed(output/'airy.npz', **arrays)
    payload = dict(airy=airy, near_mesh=near)
    source = Path('build/fullwave_imaging/default/scattering')
    from .parameters import Parameters
    from .physical_audit import sphere_rule
    from .far_field import radiation_amplitude
    p = Parameters(**json.loads((source/'report.json').read_text())['parameters'])
    with np.load(source/'paired_fields.npz') as data:
        points, weights = data['source_points_with'], data['source_polarization_weights_with']
    angular = []
    for order in (20, 32, 48):
        directions, measure = sphere_rule(order, hemisphere=True)
        amplitude = radiation_amplitude(points, weights, directions, p)
        angular.append(dict(gauss_order=order, upper_dsc_um2=float(np.sum(measure*np.sum(abs(amplitude)**2, axis=-1)))))
    payload['angular_quadrature'] = angular
    (output/'optics_and_near.json').write_text(json.dumps(payload, indent=2), encoding='utf-8')
    print(json.dumps(payload, indent=2))


if __name__ == '__main__':
    main()
