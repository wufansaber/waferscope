"""Continuous radiation of one frozen cell current, not an isolated-defect solve."""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np

from .geometry import Scene
from .contrast_current import modal_reconstruction, project_ray


def matching_scenes(left, right):
    """Allow only decimal serialization roundoff, not changed model inputs."""
    a, b = asdict(left), asdict(right)
    return all(np.isclose(value, b[key], rtol=5e-14, atol=0.) if isinstance(value, float)
               else value == b[key] for key, value in a.items())


def load_source(folder):
    folder = Path(folder)
    report = json.loads((folder/'report.json').read_text(encoding='utf-8'))
    raw = folder/'modes.npz'
    digest = hashlib.sha256(raw.read_bytes()).hexdigest()
    if report.get('raw_npz_sha256') != digest:
        raise ValueError('Raw NPZ hash does not match report: {}'.format(raw))
    with np.load(raw, allow_pickle=False) as archive:
        data = {key: archive[key] for key in archive.files}
    return data, Scene(**report['scene']), report, dict(folder=str(folder.resolve()), sha256=digest)


def cone_quadrature(na, tilt_deg, radial, azimuthal):
    if not np.isfinite(na) or not 0 <= na < 1 or not np.isfinite(tilt_deg):
        raise ValueError('Finite air NA in [0,1) and tilt required')
    if abs(tilt_deg)+np.rad2deg(np.arcsin(na)) >= 90:
        raise ValueError('This diagnostic requires the whole cone strictly above the wafer')
    if radial < 2 or azimuthal < 4:
        raise ValueError('At least 2 radial and 4 azimuthal samples required')
    nodes, weights = np.polynomial.legendre.leggauss(radial)
    lower = np.sqrt(1-na**2)
    mu = lower+(nodes+1)*(1-lower)/2
    wm = weights*(1-lower)/2
    phi = (np.arange(azimuthal)+.5)*2*np.pi/azimuthal
    angle = np.deg2rad(tilt_deg)
    axis = np.array([np.sin(angle), 0., np.cos(angle)])
    u = np.array([np.cos(angle), 0., -np.sin(angle)])
    v = np.array([0., 1., 0.])
    direction = (mu[:, None, None]*axis+np.sqrt(1-mu**2)[:, None, None]*(
        np.cos(phi)[None, :, None]*u+np.sin(phi)[None, :, None]*v))
    return direction.reshape(-1, 3), np.repeat(wm, azimuthal)*2*np.pi/azimuthal


def integrate_cones(source, scene, na, tilts, radial, azimuthal, reference_halfwidth):
    directions, weights, fields = [], [], []
    for tilt in tilts:
        rays, solid_angle = cone_quadrature(na, tilt, radial, azimuthal)
        electric = np.array([project_ray(source, scene, ray, reference_halfwidth=reference_halfwidth)
                             for ray in rays])
        directions.append(rays)
        weights.append(solid_angle)
        fields.append(electric)
    fields, weights = np.array(fields), np.array(weights)
    dsc = np.sum(abs(fields)**2, axis=-1)
    return dict(directions=np.array(directions), solid_angle_weights=weights,
                far_electric_um=fields, dsc_um2_per_sr=dsc,
                collected_cross_section_um2=np.sum(weights*dsc, axis=-1))


def lattice_diagnostic(source, scene, lattice, lattice_scene, reference_halfwidth=40,
                       na=.6, tilts=(45., -45., 0.)):
    """Change only angular sampling, with exactly the same frozen cell current."""
    valid = np.flatnonzero(lattice['propagating'])
    rays = lattice['directions'][valid]
    fields = np.array([project_ray(source, scene, ray, reference_halfwidth=reference_halfwidth)
                       for ray in rays])
    k0 = 2*np.pi/scene.wavelength_um
    area = lattice_scene.lx_um*lattice_scene.ly_um
    weights = (2*np.pi)**2/(area*k0*k0*rays[:, 2])
    dsc = np.sum(abs(fields)**2, axis=-1)
    channels, counts = [], []
    for tilt in tilts:
        angle = np.deg2rad(tilt)
        axis = np.array([np.sin(angle), 0., np.cos(angle)])
        take = (rays@axis >= np.sqrt(1-na**2)) & (na > 0)
        channels.append(float(np.sum(weights[take]*dsc[take])))
        counts.append(int(take.sum()))
    return dict(directions=rays, orders=lattice['orders'][valid], far_electric_um=fields,
                lattice_weights=weights, dsc_um2_per_sr=dsc,
                collected_cross_section_um2=np.array(channels), collected_orders=np.array(counts),
                hemisphere_cross_section_um2=np.asarray(np.sum(weights*dsc)))


def run(source_folder, reference_folder, output, levels, na=.6, tilts=(45., -45., 0.)):
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError('Use an empty output folder; existing results are never overwritten')
    levels = [tuple(level) for level in levels]
    if not levels or len(levels) != len(set(levels)) or any(len(level) != 3 or
            any(isinstance(n, bool) or not isinstance(n, int) for n in level) or
            level[0] < 2 or level[1] < 4 or level[2] < 1 for level in levels):
        raise ValueError('Distinct integer (radial>=2, azimuthal>=4, reference>=1) levels required')
    for tilt in tilts:
        cone_quadrature(na, tilt, levels[0][0], levels[0][1])
    source, scene, info, provenance = load_source(source_folder)
    reference, ref_scene, ref_info, ref_provenance = load_source(reference_folder)
    if not matching_scenes(scene, ref_scene) or not info.get('particle') or ref_info.get('particle'):
        raise ValueError('Matching particle/reference scenes are required')
    if 'source_electric_fourier' not in source or not info.get('current_z_order'):
        raise ValueError('An explicit analytic contrast-current export is required')
    if scene.kind == 'grating' and scene.grating_method != 'analytic_li':
        raise ValueError('Matching analytic Li background required')
    output.mkdir(parents=True, exist_ok=True)
    reconstructed, gate = modal_reconstruction(source, scene, reference)
    np.savez_compressed(output/'modal_gate.npz', reconstructed_electric=reconstructed,
                        direct_electric=source['electric_reflected_z0']-reference['electric_reflected_z0'],
                        directions=source['directions'], orders=source['orders'])
    (output/'gate.json').write_text(json.dumps(gate, indent=2), encoding='utf-8')
    if (not gate['comparable'] or gate['relative_power_weighted_complex_field_change'] is None or
            not np.isfinite(gate['relative_power_weighted_complex_field_change']) or
            gate['relative_power_weighted_complex_field_change'] > 1e-7 or
            gate['relative_total_power_difference'] is None or
            not np.isfinite(gate['relative_total_power_difference']) or
            abs(gate['relative_total_power_difference']) > 2e-7):
        raise RuntimeError('Same-basis absolute complex-field gate failed; no cone results emitted')
    rows, previous = [], None
    for radial, azimuthal, reference_halfwidth in levels:
        result = integrate_cones(source, scene, na, tilts, radial, azimuthal, reference_halfwidth)
        filename = 'cone_{}_{}_ref{}.npz'.format(radial, azimuthal, reference_halfwidth)
        np.savez_compressed(output/filename, **result)
        cross = result['collected_cross_section_um2']
        row = dict(radial=radial, azimuthal=azimuthal, reference_halfwidth=reference_halfwidth,
                   raw=filename, raw_sha256=hashlib.sha256((output/filename).read_bytes()).hexdigest(),
                   cross_section_um2=cross.tolist())
        if previous is not None:
            row['relative_change_from_previous'] = [float(a/b-1) if b else None
                                                    for a, b in zip(cross, previous)]
        previous = cross
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(scene=asdict(scene), source=provenance, reference=ref_provenance,
        same_basis_gate=gate, na=na, tilts_deg=list(tilts), levels=rows,
        status='FROZEN_PERIODIC_CURRENT_DIAGNOSTIC_NOT_ISOLATED_DEFECT',
        meaning='Radiation from the selected cell contrast current in a bare periodic background.',
        limitations=['Source current still contains periodic-copy interactions.',
                     'Finite Fourier truncation, z slicing, and angular quadrature need separate checks.',
                     'Not |E_sample|^2-|E_reference|^2; no camera SNR certification.'])
    (output/'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report
