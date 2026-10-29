"""Frozen PSL reference cases, separate from production UI and material defaults."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform

import numpy as np

from fullwave_prototype.parameters import Parameters


PAPER_DOI = 'https://doi.org/10.1364/AO.41.005405'
PAPER_MATERIALS = {
    442: (1.615, 4.76, .17, 1.466),
    532: (1.598, 4.15, .05, 1.461),
    633: (1.588, 3.88, .02, 1.457),
}


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')


def homogeneous_parameters(order, mesh_um=.06):
    return Parameters(wavelength_nm=266., diameter_nm=100., particle_n=1.59,
        particle_k=0., substrate_n=1., substrate_k=0., gap_nm=0.,
        incidence_deg=75., polarization='p', domain_radius_um=.18,
        pml_thickness_um=.16, pml_strength=1., pml_order=2,
        pml_mesh_scale=1.8, air_mesh_um=mesh_um, substrate_mesh_um=mesh_um,
        particle_mesh_um=.01, order=order, integration_bonus=8,
        static_condensation=True, direct_solver='sparsecholesky')


def mie_reference(output):
    import miepython
    from fullwave_prototype.theory_reference import bessel_mie_cross_sections
    p = homogeneous_parameters(2)
    exact = bessel_mie_cross_sections(p)
    qext, qsca, _, _ = miepython.efficiencies(p.particle_n, .1, .266)
    area = np.pi*.05**2
    result = dict(parameters=asdict(p), direct_bessel=exact,
        miepython_csca_um2=float(qsca*area), miepython_cext_um2=float(qext*area),
        scattering_relative_error=abs(qsca*area/exact['csca_um2']-1),
        extinction_relative_error=abs(qext*area/exact['cext_um2']-1),
        miepython_version=miepython.__version__, numpy_version=np.__version__,
        python_version=platform.python_version(),
        scope='Homogeneous air only. Legacy ideal PSL n=1.59 at 266 nm, not traced batch dispersion.')
    save_json(output/'mie_reference.json', result)
    print(json.dumps(result, indent=2))


def nist_mie(output):
    from pySCATMECH.scatterer import Free_Space_Scatterer
    from pySCATMECH.local import Local_BRDF_Model
    p = homogeneous_parameters(2)
    model = Free_Space_Scatterer('MieScatterer', {'lambda': .266, 'radius': .05,
        'sphere': complex(1.59), 'medium': 1.})
    mu, weights = np.polynomial.legendre.leggauss(64)
    theta = np.arccos(mu)
    dsc = np.array([np.asarray(model.DifferentialScatteringCrossSection(
        (0., 0., 1.), (np.sin(t), 0., np.cos(t))))[0, 0] for t in theta])
    sigma = float(2*np.pi*weights@dsc)
    result = dict(parameters=asdict(p), csca_um2=sigma,
        cext_um2=float(np.asarray(model.Extinction((0., 0., 1.)))[0, 0]),
        quadrature_order=64, note='Unpolarized DSC: sphere total cross section is polarization independent.')
    np.savez_compressed(output/'nist_mie_raw.npz', cosine=mu, weights=weights, dsc_um2_sr=dsc)
    mu, w = np.polynomial.legendre.leggauss(20)
    mu, w = (mu+1)/2, w/2
    mm, az = np.meshgrid(mu, np.arange(40)*np.pi/20, indexing='ij')
    rho = np.sqrt(1-mm**2)
    rays = np.stack((rho*np.cos(az), rho*np.sin(az), mm), axis=-1).reshape(-1, 3)
    bv = Local_BRDF_Model('Bobbert_Vlieger_BRDF_Model', {'lambda': .266, 'radius': .05,
        'sphere': complex(1.59), 'substrate': complex(1.), 'delta': 0.,
        'lmax': 12, 'order': -1, 'Norm_Inc_Approx': 0, 'improve': 3})
    jones = np.array([bv.JonesDSC(np.deg2rad(75.), np.arccos(ray[2]),
        np.arctan2(-ray[0], ray[1]), coords='psps') for ray in rays])
    np.savez_compressed(output/'nist_bv_air_raw.npz', directions=rays,
        weights_sr=np.repeat(w*np.pi/20, 40), jones_um=jones,
        p_dsc_um2_sr=np.sum(abs(jones[:, :, 1])**2, axis=1))
    save_json(output/'nist_mie.json', result)
    print(json.dumps(result, indent=2))


def paper_case(wavelength, diameter, lmax):
    from scatmech_ui.parameters import SimulationParameters, CoatingLayerParameters
    from scatmech_ui.surface_model import build_scattering_model
    pn, sn, sk, oxide_n = PAPER_MATERIALS[wavelength]
    p = SimulationParameters(wavelength_nm=float(wavelength),
        particle_diameter_nm=float(diameter), particle_n=pn, particle_k=0.,
        wafer_n=sn, wafer_k=sk, incidence_deg=60., roughness_enabled=False,
        coating_layers=(CoatingLayerParameters(n=oxide_n, k=0., thickness_nm=1.5),))
    model = build_scattering_model(p).model
    model.setParameters(delta=0., lmax=lmax, order=-1, Norm_Inc_Approx=0, improve=3)
    angles = np.arange(-89., 89.01, .5)
    # Positive signed viewing angle is the specular half-plane (SCATMECH phi=0).
    jones = np.array([model.JonesDSC(np.pi/3, np.deg2rad(abs(t)),
        0. if t >= 0 else np.pi, coords='psps') for t in angles])
    # Equal coherent s/p components, not the incoherent unpolarized state.
    field = jones @ (np.ones(2)/np.sqrt(2))
    dsc = np.sum(abs(field)**2, axis=1)
    return p, angles, jones, field, dsc


def literature(output):
    rows = []
    for wavelength, diameter in ((442, 101), (442, 155), (442, 218), (532, 155), (633, 155)):
        p, angles, j12, _, d12 = paper_case(wavelength, diameter, 12)
        _, _, j16, field, d16 = paper_case(wavelength, diameter, 16)
        name = 'paper_{}nm_d{}nm'.format(wavelength, diameter)
        raw = output/(name+'.npz')
        np.savez_compressed(raw, theta_signed_deg=angles, jones_l12_um=j12,
            jones_l16_um=j16, field_equal_sp_um=field,
            dsc_l12_um2_sr=d12, dsc_l16_um2_sr=d16)
        np.savetxt(output/(name+'.csv'), np.column_stack((angles, d16)), delimiter=',',
            header='view_angle_deg,dsc_um2_sr', comments='')
        simulation = asdict(p)
        del simulation['polarization']  # Jones excitation below overrides the UI's s/p selector.
        rows.append(dict(name=name, simulation=simulation, gap_nm=0.,
            incident_polarization='coherent equal s/p; independent Jones override, not a UI preset',
            incident_jones_sp=[float(1/np.sqrt(2))]*2,
            max_dsc_change_l12_l16_relative=float(np.max(abs(d12/d16-1))),
            jones_change_l12_l16_relative_l2=float(np.linalg.norm(j12-j16)/np.linalg.norm(j16)),
            peak_angle_deg=float(angles[np.argmax(d16)]), peak_dsc_um2_sr=float(d16.max()),
            raw_sha256=hashlib.sha256(raw.read_bytes()).hexdigest(),
            experimental_comparison_status='not_digitized_no_pointwise_error_available'))
    result = dict(doi=PAPER_DOI, cases=rows,
        note='Matched literature inputs, NOT a measured-data comparison. No fitted scale or amplitude normalization. '
             'SCATMECH and the publication share the BV implementation lineage; not an independent solver.',
        viewing_convention='Positive theta: specular half-plane phi=0. Negative theta: phi=pi. '
            'SCATMECH psps basis used on both half-planes; no published eta/Pc comparison asserted.',
        experimental_exclusions='Paper notes specular contamination near +60 deg, obscuration near -60 deg, '
            'and uncertain absolute density for the 101 nm sample.')
    save_json(output/'literature_inputs_and_runs.json', result)
    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', required=True, choices=('mie', 'nist-mie', 'literature', 'fem'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--order', type=int, choices=(2, 3), default=2)
    parser.add_argument('--mesh-um', type=float, default=.06)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.stage == 'fem':
        from fullwave_prototype.theory_audit import sphere_check
        sphere_check(homogeneous_parameters(args.order, args.mesh_um), args.output)
    else:
        {'mie': mie_reference, 'nist-mie': nist_mie, 'literature': literature}[args.stage](args.output)


if __name__ == '__main__':
    main()
