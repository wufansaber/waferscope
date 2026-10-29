"""Independent per-order comparison with a saved SCATMECH RCW benchmark."""
import hashlib
import json
from pathlib import Path

import numpy as np


def nist_comparison(scene, data, oracle_path):
    oracle_path = Path(oracle_path)
    metadata = json.loads(oracle_path.with_suffix('.json').read_text(encoding='utf-8'))
    mapping = dict(wavelength_um='wavelength_um', period_um='pitch_um',
        ridge_width_um='ridge_width_um', depth_um='depth_um', substrate_n='substrate_n',
        substrate_k='substrate_k', incidence_deg='incidence_deg')
    for key, target in mapping.items():
        if metadata['inputs'][key] != getattr(scene, target):
            raise ValueError('NIST and S4 inputs differ: '+key)
    if metadata['polarization'] != 'p' or scene.kind != 'grating':
        raise ValueError('This oracle only validates a bare grating with p illumination')
    with np.load(oracle_path, allow_pickle=False) as saved:
        oracle = {key: saved[key] for key in saved.files}
    lookup = {tuple(g): i for i, g in enumerate(data['orders'])}
    per_order, numerator, denominator = [], 0., 0.
    for order, direction, amplitude, expected in zip(oracle['orders'], oracle['direction'],
                                                   oracle['amplitude'], oracle['efficiency']):
        if expected <= 1e-12:
            continue
        index = lookup.get((0, int(order)*scene.periods_y))
        if index is None or not data['propagating'][index]:
            raise ValueError('S4 basis omits a propagating benchmark order')
        # NIST transverse +x is application +KY; its p vector rotates with the ray.
        p = np.array([0., -direction[2], direction[0]])
        expected_e = amplitude*p
        actual_e = data['electric_reflected_z0'][index]
        actual_p = np.dot(actual_e, p)
        weight = direction[2]/np.cos(np.deg2rad(scene.incidence_deg))
        numerator += weight*np.sum(abs(actual_e-expected_e)**2)
        denominator += expected
        actual = float(data['reflectance_by_order'][index])
        np.testing.assert_allclose(data['directions'][index], [0., direction[0], direction[2]], atol=1e-12)
        per_order.append(dict(order=int(order), nist_efficiency=float(expected),
            s4_efficiency=actual, absolute_power_error=actual-float(expected),
            relative_power_error=actual/expected-1,
            nist_p_amplitude=[float(amplitude.real), float(amplitude.imag)],
            s4_p_amplitude=[float(actual_p.real), float(actual_p.imag)],
            complex_amplitude_absolute_error=float(abs(actual_p-amplitude)),
            phase_error_deg=float(np.rad2deg(np.angle(actual_p*np.conj(amplitude))))))
    return dict(oracle=str(oracle_path), oracle_sha256=hashlib.sha256(oracle_path.read_bytes()).hexdigest(),
        oracle_inputs=metadata, per_order=per_order,
        relative_total_reflectance_error=float(data['reflectance_by_order'].sum()/oracle['efficiency'].sum()-1),
        power_weighted_relative_complex_field_error=float(np.sqrt(numerator/denominator)),
        note='No amplitude fit, phase fit, or symmetry enforcement. Bare periodic grating only.')
